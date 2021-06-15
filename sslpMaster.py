import numpy as np
import ast
from gurobipy import *
import master
import time

class SSLPinst(master.StochIPinst):
	def __init__(self,instance):
		self.Nscen = None
		self.Nfv = None
		self.Nn = None
		self.instance = instance
		self.name = instance
		self.cVec = None
		self.qMtr = None
		self.q0Vec = None
		self.dMtr = None
		self.Nu = None
		self.hMtr = None
		self.pVec = None
		self.PrimalMaster = Model('Master')
		self.theta = {}
		self.x = {}
		# Benders cuts lists
		self.cutlist = {}
		# All cuts coef. lists
		self.coeflist = {}
		self.thetaCutList = {}
		self.Ncuts = 0
		self.Nsubs = 0
		# Subproblems
		self.scensub = {}
		self.scenrc = {}

	def _readVec(self, vec, nType=np.float32):
		return np.asarray(ast.literal_eval(vec[:-1]), dtype = nType)

	def _readMtr(self, mtr):
		lists_string = mtr.strip('][\n').split('; ')
		lists_array = []
		for line in lists_string:
			lists_array.append(line.split(' '))
		return np.asarray(lists_array, dtype = np.float32)

	def Initialize(self):
		objCoef = {}
		for i in range(self.Nfv):
			objCoef[i] = self.cVec[i]
		regCoefy = 1.0
		t0 = time.time()
		for s in range(self.Nscen):
			self.scensub[s] = Model()
			self.scensub[s].setParam('OutputFlag', False)
			self.scensub[s].modelSense = GRB.MINIMIZE
			x_scen = {}
			for j in range(self.Nfv):
				x_scen[j] = self.scensub[s].addVar(vtype=GRB.BINARY,obj=objCoef[j],name="x"+str(j))
			y = {}
			for i in range(self.Nn):
				for j in range(self.Nfv):
					y[i,j] = self.scensub[s].addVar(vtype=GRB.BINARY,obj=-regCoefy*self.qMtr[i,j],name="y_"+str(i)+","+str(j))
			y0 = {}
			for j in range(self.Nfv):
				y0[j] = self.scensub[s].addVar(lb=0.0,obj=regCoefy*self.q0Vec[j],name="y0_"+str(j))
			for j in range(self.Nfv):
				self.scensub[s].addConstr(quicksum(self.dMtr[i,j]*y[i,j] for i in range(self.Nn))-y0[j] <= self.Nu*x_scen[j])
			for i in range(self.Nn):
				self.scensub[s].addConstr(quicksum(y[i,j] for j in range(self.Nfv)) == self.hMtr[i,s])
			self.scensub[s].update()
		print("Build scenario subproblems: "+str(time.time()-t0))
		t0 = time.time()
		for s in range(self.Nscen):
			self.scenrc[s] = Model()
			self.scenrc[s].setParam('OutputFlag', False)
			self.scenrc[s].modelSense = GRB.MINIMIZE
			for i in range(self.Nn):
				for j in range(self.Nfv):
					y[i,j] = self.scenrc[s].addVar(vtype=GRB.BINARY,obj=-self.qMtr[i,j],name="y_"+str(i)+","+str(j))
			y0 = {}
			for j in range(self.Nfv):
				y0[j] = self.scenrc[s].addVar(lb=0.0,obj=self.q0Vec[j],name="y0_"+str(j))
			for j in range(self.Nfv):
				self.scenrc[s].addConstr(quicksum(self.dMtr[i,j]*y[i,j] for i in range(self.Nn))-y0[j] <= 0.0,name="xconstr"+str(j))
			for i in range(self.Nn):
				self.scenrc[s].addConstr(quicksum(y[i,j] for j in range(self.Nfv)) == self.hMtr[i,s])
			self.scenrc[s].update()
		print("Build recourse subproblems: "+str(time.time()-t0))

	def FreeMemory(self):
		for s in range(self.Nscen):
			del self.scensub[s]
			del self.scenrc[s]
		del self.PrimalMaster


	def SolveExtensive(self):
		Extensive = Model('Extensive')
		x = {}
		y = {}
		y0 = {}
		for j in range(self.Nfv):
			x[j] = Extensive.addVar(vtype=GRB.BINARY,obj=self.cVec[j])
		for s in range(self.Nscen):
			for i in range(self.Nn):
				for j in range(self.Nfv):
					y[s,i,j] = Extensive.addVar(vtype=GRB.BINARY,obj=-self.pVec[s]*self.qMtr[i,j])
		for s in range(self.Nscen):
			for j in range(self.Nfv):
				y0[s,j] = Extensive.addVar(lb=0.0,obj=self.pVec[s]*self.q0Vec[j])
		for s in range(self.Nscen):
			for j in range(self.Nfv):
				Extensive.addConstr(quicksum(self.dMtr[i,j]*y[s,i,j] for i in range(self.Nn))-y0[s,j] <= self.Nu*x[j])
		for s in range(self.Nscen):
			for i in range(self.Nn):
				Extensive.addConstr(quicksum(y[s,i,j] for j in range(self.Nfv)) == self.hMtr[i,s])
		Extensive.modelSense = GRB.MINIMIZE
		Extensive.update()
		t0 = time.time()
		Extensive.optimize()
		print('Solve Extensive Form: '+str(self.name))
		print('Optimal Obj. Value: '+str(Extensive.objval)+', Total Time: '+str(time.time()-t0))
		x_value = {}
		for i in range(self.Nfv):
		    x_value[i] = x[i].x
		print(x_value)
		theta_value = {}
		for s in range(self.Nscen):
			theta_value[s] = sum(-self.qMtr[i,j]*y[s,i,j].x for i in range(self.Nn) for j in range(self.Nfv))+sum(self.q0Vec[j]*y0[s,j].x for j in range(self.Nfv))
		print(theta_value)


	def BuildPrimal(self):
		for s in range(self.Nscen):
			self.theta[s] = self.PrimalMaster.addVar(lb=-GRB.INFINITY)
		for i in range(self.Nfv):
			self.x[i] = self.PrimalMaster.addVar(lb=0.0,ub=1.0)
		self.PrimalMaster.setObjective(quicksum(self.pVec[s]*self.theta[s] for s in range(self.Nscen))+quicksum(self.cVec[j]*self.x[j] for j in range(self.Nfv)))
		self.PrimalMaster.setParam( 'OutputFlag', False )
		self.PrimalMaster.modelSense = GRB.MINIMIZE
		self.PrimalMaster.update()
		self.addBenders(init=True)

	def readData(self):
		f = open('Instances/sslp/'+self.instance+'.txt', 'r')
		x = f.readlines()
		f.close()
		fline = self._readVec(x[0], np.int32)
		self.Nfv = fline[0]
		self.Nn = fline[1]
		self.Nscen = fline[2]
		# Read cVec
		self.cVec = self._readVec(x[1], np.float32)
		# Read qMtr
		self.qMtr = self._readMtr(x[2])
		# Read q0Vec
		self.q0Vec = self._readVec(x[3], np.float32)
		# Read dMtr
		self.dMtr = self._readMtr(x[4])
		# Read Nu
		self.Nu = float(x[5])
		#Read hMtr
		self.hMtr = self._readMtr(x[6])
		#Read pVec
		self.pVec = self._readVec(x[7], np.float32)
		for s in range(self.Nscen):
			self.cutlist[s] = []
			self.coeflist[s] = []
			self.thetaCutList[s] = []
		self.Initialize()
		
	def SolveBendersSub(self,scen_id,x_input):
		benderssub = self.scenrc[scen_id].relax()
		for j in range(self.Nfv):
			benderssub.getConstrByName("xconstr"+str(j)).RHS = self.Nu*x_input[j]
		benderssub.update()
		benderssub.optimize()


		subg = {}
		for j in range(self.Nfv):
			subg[j] = self.Nu*benderssub.getConstrByName("xconstr"+str(j)).pi
		const = benderssub.objval-sum(subg[j]*x_input[j] for j in range(self.Nfv))
		return benderssub.objval,const,subg

	def SolveScenSub(self,scen_id,objCoef,regCoefy,tlimit,MIPgap='nan'):
		self.scensub[scen_id].params.TimeLimit = tlimit
		if MIPgap != 'nan':
			self.scensub[scen_id].params.MIPGap = MIPgap
		for j in range(self.Nfv):
			self.scensub[scen_id].getVarByName("x"+str(j)).obj = objCoef[j]
		for i in range(self.Nn):
			for j in range(self.Nfv):
				self.scensub[scen_id].getVarByName("y_"+str(i)+","+str(j)).obj = -regCoefy*self.qMtr[i,j]
		for j in range(self.Nfv):
			self.scensub[scen_id].getVarByName("y0_"+str(j)).obj = regCoefy*self.q0Vec[j]
		self.scensub[scen_id].update()
		self.scensub[scen_id].optimize()
		if self.scensub[scen_id].status != 2:
			print('IP Status: '+str(self.scensub[scen_id].status))
		ObjReturn = self.scensub[scen_id].objval
		BestBound = self.scensub[scen_id].ObjBound
		x_value = {}
		for j in range(self.Nfv):
			x_value[j] = self.scensub[scen_id].getVarByName("x"+str(j)).x
		yObjV = -sum(self.qMtr[i,j]*self.scensub[scen_id].getVarByName("y_"+str(i)+","+str(j)).x for i in range(self.Nn) for j in range(self.Nfv))+\
			sum(self.q0Vec[j]*self.scensub[scen_id].getVarByName("y0_"+str(j)).x for j in range(self.Nfv))

		self.Nsubs += 1
		SubOpt = []
		for k in range(self.scensub[scen_id].SolCount-1):
			self.scensub[scen_id].setParam('SolutionNumber', k+1)
			soln = {}
			modelSoln = self.scensub[scen_id].Xn
			for j in range(self.Nfv):
				soln[j] = modelSoln[j]
			soln[self.Nfv]=-sum(self.qMtr[i,j]*modelSoln[self.Nfv+i*self.Nfv+j] for i in range(self.Nn) for j in range(self.Nfv))+\
				sum(self.q0Vec[j]*modelSoln[self.Nfv+self.Nn*self.Nfv+j] for j in range(self.Nfv))
			SubOpt.append(soln)

		if regCoefy < 1e-4:
			yObjV = self.SolveRecourse(scen_id,x_value)
			ObjReturn = sum(objCoef[j]*x_value[j] for j in range(self.Nfv))+regCoefy*yObjV

		return (ObjReturn,x_value,yObjV,SubOpt,BestBound)

	def SolveRecourse(self,scen_id,x_input):
		for j in range(self.Nfv):
			self.scenrc[scen_id].getConstrByName("xconstr"+str(j)).RHS = self.Nu*x_input[j]
		self.scenrc[scen_id].update()
		self.scenrc[scen_id].optimize()
		if self.scenrc[scen_id].status != 2:
			print('IP Status: '+str(self.scenrc[scen_id].status))
		yObjV = self.scenrc[scen_id].objval
		return yObjV

	def BendersBnC(self,timeLimit=60*60,tol=1e-6,cutTime=0,methodName='BnC'):
		t0 = time.time()
		global tBenders
		global tRecourse
		global MIPSOLcount
		tBenders = 0.0
		tRecourse = 0.0
		MIPSOLcount = 0
		ScenLB = {}
		for s in range(self.Nscen):
			ScenLB[s] = -sum(max(self.qMtr[i,j] for j in range(self.Nfv))*self.hMtr[i,s] for i in range(self.Nn))
		for i in range(self.Nfv):
			self.x[i].vtype = GRB.BINARY
		self.PrimalMaster.setParam( 'OutputFlag', True )
		self.PrimalMaster.Params.LogFile="Results/"+str(self.name)+"_"+methodName+".log"
		self.PrimalMaster.update()
		def cb(model, where):
			if where == GRB.Callback.MIPSOL:
				global tBenders
				global tRecourse
				global MIPSOLcount
				MIPSOLcount = MIPSOLcount+1
				x_value = {}
				theta_value = {}
				J1 = []
				J0 = []
				for i in range(self.Nfv):
					x_value[i] = model.cbGetSolution(self.x[i])
					if x_value[i] > 0.5:
						J1.append(i)
						x_value[i] = 1.0
					else:
						J0.append(i)
						x_value[i] = 0.0
				for s in range(self.Nscen):
					theta_value[s] = model.cbGetSolution(self.theta[s])
				for s in range(self.Nscen):
					tRecourseStart = time.time()
					ObjV = self.SolveRecourse(scen_id=s,x_input=x_value)
					tRecourse = tRecourse+(time.time()-tRecourseStart)
					if theta_value[s] < ObjV-tol*(abs(theta_value[s])+1):
						model.cbLazy(self.theta[s] >= ObjV-(ObjV-ScenLB[s])*(quicksum((1-self.x[j]) for j in J1)+quicksum(self.x[j] for j in J0)))
					tBendersStart = time.time()
					ObjV,const,subg = self.SolveBendersSub(scen_id=s,x_input=x_value)
					tBenders = tBenders+(time.time()-tBendersStart)
					if theta_value[s] < ObjV-tol*(abs(theta_value[s])+1):
						model.cbLazy(self.theta[s] >= const+quicksum(subg[i]*self.x[i] for i in range(self.Nfv)))
		self.PrimalMaster.Params.Heuristics = 0.0
		self.PrimalMaster.Params.lazyConstraints = 1
		self.PrimalMaster.Params.timeLimit = timeLimit-(time.time()-t0)
		print("Time_Limit: "+str(timeLimit-(time.time()-t0)))
		self.PrimalMaster.optimize(cb)
		
		print("Optimal Value: "+str(self.PrimalMaster.objval)+", time: "+str(time.time()-t0))
		wrtStr = "ObjVal: "+str(self.PrimalMaster.objval)+"\tObjBnd: "+str(self.PrimalMaster.ObjBound)+"\tBnCTime: "+str(time.time()-t0)+"\tNo. of MIP solns: "+str(MIPSOLcount)+"\tBendersTime: "+str(tBenders)+"\tRecourseTime: "+str(tRecourse)+"\tTotalTime: "+str(cutTime+time.time()-t0)+"\n"
		f = open("Results/"+str(self.name)+"_"+methodName+".txt","a")
		f.write(wrtStr)
		f.close