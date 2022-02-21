import numpy as np
import ast
from gurobipy import *
import master
import time

class SNIPinst(master.StochIPinst):
	def __init__(self,instance,R,snipNo):
		self.instance = instance
		self.name = str(instance)+'_'+str(R)+'_'+str(snipNo)
		self.R = R
		self.snipNo = snipNo
		self.Nscen = None
		self.Nfv = None
		self.cVec = None
		self.pVec = None
		self.A = []
		self.AD = []
		self.ALL = []
		self.N = []
		self.SCEN = []
		self.PSI = []
		self.ind = {}
		self.n = None
		self.na = None
		self.a = None
		self.ad = None
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
		self.scenrc = None

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
			for i in range(self.Nfv):
				x_scen[i] = self.scensub[s].addVar(vtype=GRB.BINARY,obj=objCoef[i],name="x"+str(i))
			ypi = {}
			for i in range(self.n):
				ypi[i] = self.scensub[s].addVar(lb=0.0,ub=1.0,obj=0.0,name="y"+str(i))
			ypi[self.ind[self.SCEN[s,0]]].obj = regCoefy
			self.scensub[s].addConstr(quicksum(x_scen[i] for i in range(self.Nfv)) <= self.R)
			for i in range(self.a):
				self.scensub[s].addConstr(ypi[self.ind[self.A[i][0]]] - self.A[i][2]*ypi[self.ind[self.A[i][1]]] >= 0)
			for i in range(self.ad):
				self.scensub[s].addConstr(ypi[self.ind[self.AD[i][0]]] - self.AD[i][2]*ypi[self.ind[self.AD[i][1]]] >= -(self.AD[i][2] - self.AD[i][3])*self.PSI[s,self.ind[self.AD[i][1]]]*x_scen[i])
				self.scensub[s].addConstr(ypi[self.ind[self.AD[i][0]]] - self.AD[i][3]*ypi[self.ind[self.AD[i][1]]] >= 0)
			self.scensub[s].addConstr(ypi[self.ind[self.SCEN[s,1]]] == 1)
			self.scensub[s].setParam( 'OutputFlag', False )
			self.scensub[s].modelSense = GRB.MINIMIZE
			self.scensub[s].update()
		print("Build scenario subproblems: "+str(time.time()-t0))
		t0 = time.time()
		self.scenrc = Model()
		self.scenrc.setParam('OutputFlag', False)
		self.scenrc.modelSense = GRB.MINIMIZE
		ypi = {}
		subgConstr = {}
		for i in range(self.n):
			ypi[i] = self.scenrc.addVar(lb=0.0,ub=1.0,obj=0.0,name="y"+str(i))
		for i in range(self.a):
			self.scenrc.addConstr(ypi[self.ind[self.A[i][0]]] - self.A[i][2]*ypi[self.ind[self.A[i][1]]] >= 0)
		for i in range(self.ad):
			subgConstr[i] = self.scenrc.addConstr(ypi[self.ind[self.AD[i][0]]] - self.AD[i][2]*ypi[self.ind[self.AD[i][1]]] >= 0,name="subgCon"+str(i))
			self.scenrc.addConstr(ypi[self.ind[self.AD[i][0]]] - self.AD[i][3]*ypi[self.ind[self.AD[i][1]]] >= 0)
		self.scenrc.update()
		print("Build recourse subproblems: "+str(time.time()-t0))

	def FreeMemory(self):
		for s in range(self.Nscen):
			del self.scensub[s]
		del self.scenrc
		del self.PrimalMaster

	def SolveExtensive(self):
		Extensive = Model('Extensive')
		Extensive.modelSense = GRB.MINIMIZE
		x = {}
		ypi = {}
		for i in range(self.Nfv):
			x[i] = Extensive.addVar(vtype=GRB.BINARY)
		for s in range(self.Nscen):
			for i in range(self.n):
				ypi[s,i] = Extensive.addVar(lb=0.0)
		Extensive.addConstr(quicksum(x[i] for i in range(self.Nfv)) <= self.R)
		for s in range(self.Nscen):
			for i in range(self.a):
				Extensive.addConstr(ypi[s,self.ind[self.A[i][0]]] - self.A[i][2]*ypi[s,self.ind[self.A[i][1]]] >= 0)
			for i in range(self.ad):
				Extensive.addConstr(ypi[s,self.ind[self.AD[i][0]]] - self.AD[i][2]*ypi[s,self.ind[self.AD[i][1]]] >= -(self.AD[i][2] - self.AD[i][3])*self.PSI[s,self.ind[self.AD[i][1]]]*x[i])
				Extensive.addConstr(ypi[s,self.ind[self.AD[i][0]]] - self.AD[i][3]*ypi[s,self.ind[self.AD[i][1]]] >= 0)
			Extensive.addConstr(ypi[s,self.ind[self.SCEN[s,1]]] == 1)
		Extensive.setObjective(quicksum(self.pVec[s]*ypi[s,self.ind[self.SCEN[s,0]]] for s in range(self.Nscen)))
		Extensive.update()
		t0 = time.time()
		Extensive.optimize()
		print('Solve Extensive Form: '+str(self.name))
		print('Optimal Obj. Value: '+str(Extensive.objval)+', Total Time: '+str(time.time()-t0))



	def BuildPrimal(self):
		for s in range(self.Nscen):
			self.theta[s] = self.PrimalMaster.addVar(lb=-GRB.INFINITY,name="theta"+str(s))
		for i in range(self.Nfv):
			self.x[i] = self.PrimalMaster.addVar(lb=0.0,ub=1.0,name="x"+str(i))
		self.PrimalMaster.addConstr(quicksum(self.x[i] for i in range(self.Nfv)) <= self.R)
		self.PrimalMaster.setObjective(quicksum(self.pVec[s]*self.theta[s] for s in range(self.Nscen))+quicksum(self.cVec[j]*self.x[j] for j in range(self.Nfv)))
		self.PrimalMaster.setParam( 'OutputFlag', False )
		self.PrimalMaster.modelSense = GRB.MINIMIZE
		self.PrimalMaster.update()
		self.addBenders(init=True)

	def readData(self):
		instanceName = self.instance+'.txt'
		arcFile = "Instances/SNIP/nonint_"+instanceName
		intArcFile = "Instances/SNIP/int_"+instanceName
		f = open(arcFile)
		x = f.readlines()
		f.close()
		for line in x:
			if x.index(line)%2 == 0:
				z = line.strip('\n').split('\t\t')
				z = [int(z[0]),int(z[1]),float(z[2]),float(z[2])]
				self.N.append(int(z[0]))
				self.N.append(int(z[1]))
				# Unsorted
				self.A.append(z)

		f = open(intArcFile)
		x = f.readlines()
		f.close()
		for line in x:
			if x.index(line)%2 == 0:
				z = line.strip('\n').split('\t\t')
				if self.snipNo == 1:
					z = [int(z[0]),int(z[1]),float(z[2]),float(z[3])]
				elif self.snipNo == 2:
					z = [int(z[0]),int(z[1]),float(z[2]),0.5*float(z[2])]
				elif self.snipNo == 3:
					z = [int(z[0]),int(z[1]),float(z[2]),0.1*float(z[2])]
				elif self.snipNo == 4:
					z = [int(z[0]),int(z[1]),float(z[2]),0.0]
				# Unsorted
				self.N.append(int(z[0]))
				self.N.append(int(z[1]))
				self.AD.append(z)

		self.N = np.sort(np.unique(np.array(self.N)))
		self.ALL = self.A+self.AD

		self.a = np.size(self.A,0)
		self.ad = np.size(self.AD,0)
		self.Nfv = self.ad
		self.cVec = [0.0]*self.Nfv
		self.na = np.size(self.ALL,0)
		self.n = np.size(self.N,0)

		f = open("Instances/SNIP/Scenarios.txt")
		x = f.readlines()
		f.close()
		for line in x:
			z = line.strip('\n').split('\t')
			z = [int(z[0]),int(z[1]),float(z[2])]
			self.SCEN.append(z)
		self.SCEN = np.array(self.SCEN)
		self.pVec = self.SCEN[:,2]

		f = open("Instances/SNIP/psi_reformat.txt")
		x = f.readlines()
		f.close()

		for string in x:
			z = string.strip('\n').split(' ')
			z = [float(i) for i in z]
			self.PSI.append(z)
		self.PSI = np.array(self.PSI)

		for i in range(self.n):
			self.ind[self.N[i]] = i

		ns = np.size(self.SCEN,0)
		self.Nscen = ns
		for s in range(self.Nscen):
			self.cutlist[s] = []
			self.coeflist[s] = []
			self.thetaCutList[s] = []

		self.Initialize()


	def SolveBendersSub(self,scen_id,x_input):
		for i in range(self.ad):
			self.scenrc.getConstrByName("subgCon"+str(i)).RHS = -(self.AD[i][2] - self.AD[i][3])*self.PSI[scen_id,self.ind[self.AD[i][1]]]*x_input[i]
		self.scenrc.getVarByName("y"+str(self.ind[self.SCEN[scen_id,1]])).lb = 1.0
		self.scenrc.getVarByName("y"+str(self.ind[self.SCEN[scen_id,0]])).obj = 1.0
		self.scenrc.update()
		self.scenrc.optimize()
		ObjValue = self.scenrc.objval
		subg = {}
		for i in range(self.ad):
			subg[i] = -(self.AD[i][2] - self.AD[i][3])*self.PSI[scen_id,self.ind[self.AD[i][1]]]*self.scenrc.getConstrByName("subgCon"+str(i)).pi
		const = self.scenrc.objval - sum(subg[i]*x_input[i] for i in range(self.Nfv))
		self.scenrc.getVarByName("y"+str(self.ind[self.SCEN[scen_id,1]])).lb = 0.0
		self.scenrc.getVarByName("y"+str(self.ind[self.SCEN[scen_id,0]])).obj = 0.0
		self.scenrc.update()
		Benderscoef = []
		for i in range(self.Nfv):
			Benderscoef.append(subg[i])
		Benderscoef.append(const)
		Benderscoef = np.array([Benderscoef])
		self.BendersCuts[scen_id] = np.append(self.BendersCuts[scen_id], Benderscoef, axis=0)
		return ObjValue,const,subg

	def SolveScenSub(self,scen_id,objCoef,regCoefy,tlimit,MIPgap='nan'):
		self.scensub[scen_id].params.TimeLimit = tlimit
		if MIPgap != 'nan':
			self.scensub[scen_id].params.MIPGap = MIPgap
		for i in range(self.Nfv):
			self.scensub[scen_id].getVarByName("x"+str(i)).obj = objCoef[i]
		self.scensub[scen_id].getVarByName("y"+str(self.ind[self.SCEN[scen_id,0]])).obj = regCoefy
		self.scensub[scen_id].update()
		self.scensub[scen_id].optimize()
		if self.scensub[scen_id].status != 2:
			print('IP Status: '+str(self.scensub[scen_id].status))
		ObjReturn = self.scensub[scen_id].objval
		BestBound = self.scensub[scen_id].ObjBound
		x_value = {}
		for i in range(self.Nfv):
			x_value[i] = self.scensub[scen_id].getVarByName("x"+str(i)).x
		yObjV = self.scensub[scen_id].getVarByName("y"+str(self.ind[self.SCEN[scen_id,0]])).x
		self.Nsubs += 1
		SubOpt = []
		for k in range(self.scensub[scen_id].SolCount-1):
			self.scensub[scen_id].setParam('SolutionNumber', k+1)
			soln = {}
			modelSoln = self.scensub[scen_id].Xn
			for i in range(self.Nfv):
				soln[i] = modelSoln[i]
			soln[self.Nfv] = modelSoln[self.Nfv+self.ind[self.SCEN[scen_id,0]]]
			SubOpt.append(soln)

		if regCoefy < 1e-4:
			yObjV,const,subg = self.SolveBendersSub(scen_id,x_value)
			ObjReturn = sum(objCoef[j]*x_value[j] for j in range(self.Nfv))+regCoefy*yObjV
		
		return (ObjReturn,x_value,yObjV,SubOpt,BestBound)

	def BendersBnC(self,timeLimit=60*60,tol=1e-6,cutTime=0,methodName='BnC',doubleRun=False):
		t0 = time.time()
		global tBenders
		global MIPSOLcount
		tBenders = 0.0
		MIPSOLcount = 0
		for i in range(self.Nfv):
			self.x[i].vtype = GRB.BINARY
		self.PrimalMaster.setParam( 'OutputFlag', True )
		self.PrimalMaster.Params.LogFile="Results/"+str(self.name)+"_"+methodName+".log"
		self.PrimalMaster.update()

		def cb(model, where):
			if where == GRB.Callback.MIPSOL:
				global tBenders
				global MIPSOLcount
				MIPSOLcount = MIPSOLcount+1
				x_value = {}
				theta_value = {}
				for i in range(self.Nfv):
					#x_value[i] = model.cbGetSolution(self.x[i])
					x_value[i] = model.cbGetSolution(model.getVarByName("x"+str(i)))
				for s in range(self.Nscen):
					#theta_value[s] = model.cbGetSolution(self.theta[s])
					theta_value[s] = model.cbGetSolution(model.getVarByName("theta"+str(s)))
				for s in range(self.Nscen):
					tBendersStart = time.time()
					ObjV,const,subg = self.SolveBendersSub(scen_id=s,x_input=x_value)
					tBenders = tBenders+(time.time()-tBendersStart)
					if theta_value[s] < ObjV-tol*(abs(theta_value[s])+1):
						#model.cbLazy(self.theta[s] >= const+quicksum(subg[i]*self.x[i] for i in range(self.Nfv)))
						model.cbLazy(model.getVarByName("theta"+str(s)) >= const+quicksum(subg[i]*model.getVarByName("x"+str(i)) for i in range(self.Nfv)))
		self.PrimalMaster.Params.Heuristics = 0.0
		self.PrimalMaster.update()
		self.PrimalMaster.Params.lazyConstraints = 1
		self.PrimalMaster.Params.timeLimit = timeLimit-(time.time()-t0)
		print("Time_Limit: "+str(timeLimit-(time.time()-t0)))
		if doubleRun == True:
			MasterCopy = self.PrimalMaster.copy()
			MasterCopy.Params.LogFile="Results/"+str(self.name)+"_"+methodName+"_basic.log"
			MasterCopy.Params.Presolve = 0
			MasterCopy.Params.Cuts = 0
		TimeUsed = time.time()-t0
		t0 = time.time()
		self.PrimalMaster.optimize(cb)
		print("Optimal Value: "+str(self.PrimalMaster.objval)+", time: "+str(TimeUsed+time.time()-t0))
		wrtStr = "ObjVal: "+str(self.PrimalMaster.objval)+"\tObjBnd: "+str(self.PrimalMaster.ObjBound)+"\tBnCTime: "+str(TimeUsed+time.time()-t0)+"\tNo. of MIP solns: "+str(MIPSOLcount)+"\tBendersTime: "+str(tBenders)+"\tTotalTime: "+str(cutTime+TimeUsed+time.time()-t0)+"\n"
		f = open("Results/"+str(self.name)+"_"+methodName+".txt","a")
		f.write(wrtStr)
		f.close
		if doubleRun == True:
			t0 = time.time()
			tBenders = 0.0
			MIPSOLcount = 0
			MasterCopy.optimize(cb)
			print("Optimal Value: "+str(MasterCopy.objval)+", time: "+str(TimeUsed+time.time()-t0))
			wrtStr = "*ObjVal: "+str(MasterCopy.objval)+"\tObjBnd: "+str(MasterCopy.ObjBound)+"\tBnCTime: "+str(TimeUsed+time.time()-t0)+"\tNo. of MIP solns: "+str(MIPSOLcount)+"\tBendersTime: "+str(tBenders)+"\tTotalTime: "+str(cutTime+TimeUsed+time.time()-t0)+"\n"
			f = open("Results/"+str(self.name)+"_"+methodName+".txt","a")
			f.write(wrtStr)
			f.close
		