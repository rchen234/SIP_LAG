import master
import sslpMaster
import snipMaster
import time

# method = 1: Benders cuts
# method = 2: Strengthened Benders cuts
# method = 3: Exact separation of Lagrangian cuts
# method = 4: Restricted separation of Lagrangian cuts
# ifBnC: exact solution or not

def RunSnip(instanceName,Rbudget,snipNumber,method,BendersMethod='cutpl',pi0Coefficient=1e-1,tolerance=1e-6,gapTolCoef=5e-1,tLimit=60*60,maxNum=10,newn=True,mip=True,pplus=False,LevelParam=0.01,ifBnC=False,tTotalLimit=2*60*60,earlyTmn=False,delta0=0.01,doubleRun=False):
	ss = snipMaster.SNIPinst(instanceName,R=Rbudget,snipNo=snipNumber)
	ss.readData()
	ss.BuildPrimal()
	t0 = time.time()
	if method == 1:
		ss.addBenders(method=BendersMethod,store=True,ifBnC=ifBnC)
	elif method == 2:
		ss.PureStrBenders(BDmethod=BendersMethod,timeLimit=tLimit)
	elif method == 3:
		ss.IterativeLag(BDmethod=BendersMethod,pi0Coef=pi0Coefficient,tol=tolerance,gapTol=gapTolCoef,timeLimit=tLimit)
	elif method == 4:
		ss.IterRstrLag(BDmethod=BendersMethod,pi0Coef=pi0Coefficient,tol=tolerance,gapTol=gapTolCoef,timeLimit=tTotalLimit,maxNcuts=maxNum,newn=newn,mip=mip,poolplus=pplus,LevelParam=LevelParam,earlyTmn=earlyTmn)
	if ifBnC == True:
		tCut = time.time()-t0
		if method == 1:
			methodName = 'BendersBnC'
		if method == 4:
			methodName = 'BnC'
		ss.BendersBnC(timeLimit=tTotalLimit-tCut,cutTime=tCut,methodName=methodName,doubleRun=doubleRun)
	ss.FreeMemory()

def RunSslp(instanceName,method,BendersMethod='level',pi0Coefficient=1.0,tolerance=1e-6,gapTolCoef=5e-1,tLimit=60*60,maxNum=10,newn=True,mip=True,pplus=False,LevelParam=0.01,ifBnC=False,tTotalLimit=2*60*60,earlyTmn=False,delta0=0.01):
	ss = sslpMaster.SSLPinst(instanceName)
	ss.readData()
	ss.BuildPrimal()
	t0 = time.time()
	if method == 1:
		ss.addBenders(method=BendersMethod,store=True,ifBnC=ifBnC)
	elif method == 2:
		ss.PureStrBenders(BDmethod=BendersMethod,timeLimit=tLimit)
	elif method == 3:
		ss.IterativeLag(BDmethod=BendersMethod,pi0Coef=pi0Coefficient,tol=tolerance,gapTol=gapTolCoef,timeLimit=tLimit)
	elif method == 4:
		ss.IterRstrLag(BDmethod=BendersMethod,pi0Coef=pi0Coefficient,tol=tolerance,gapTol=gapTolCoef,timeLimit=tTotalLimit,maxNcuts=maxNum,newn=newn,mip=mip,poolplus=pplus,LevelParam=LevelParam,earlyTmn=earlyTmn)
	if ifBnC == True:
		tCut = time.time()-t0
		if method == 1:
			methodName = 'BendersBnC'
		if method == 4:
			methodName = 'BnC'
		ss.BendersBnC(timeLimit=tTotalLimit-tCut,cutTime=tCut,methodName=methodName)
	ss.FreeMemory()


if __name__ == '__main__':
	# Solve SNIP
	instance = 'instance0'
	Rb = 30
	snipNo = 3
	RunSnip(instance,Rb,snipNo,method=2)
	# Solve SSLP
	#RunSslp('sslp1_50_40_50')