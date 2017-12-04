#!/usr/bin/env python
import ROOT
import math
from math import *

def four_vector(obj):
	v = ROOT.TLorentzVector()
	v.SetPtEtaPhiE(obj.pt(),obj.eta(),obj.phi(),obj.e())
	return v

def pi():
	return 3.14159265359
	
def getpt(lep):
	return lep.pt()
	
def sorter(alist):
	alist.sort(key=getpt,reverse=True)

GeV = 1000.0

ROOT.gROOT.Macro( '$ROOTCOREDIR/scripts/load_packages.C' )
if(not ROOT.xAOD.Init().isSuccess()): print "Failed xAOD.Init()"

file1 = ROOT.TFile("hists.root","recreate")

#Block defining histograms

ll_mass = ROOT.TH1F("dilep_mass","M(ll)", 200,0,1000*GeV)
ll_pt = ROOT.TH1F("dilep_pt","p_{T}(ll)",200,0,1000*GeV)
ll_eta = ROOT.TH1F("dilep_eta","#eta(ll)",25,-5,5)
vv_mass = ROOT.TH1F("vv_mass","M(#nu#nu)", 200,0,1000*GeV)
vv_pt = ROOT.TH1F("vv_pt","p_{T}(#nu#nu)",200,0,1000*GeV)
vv_eta = ROOT.TH1F("vv_eta","#eta(#nu#nu)",25,-5,5)

Z_mass = ROOT.TH1F("Z_mass","M_{Z}", 200,0,1000*GeV)
M_ZZ = ROOT.TH1F("M(ZZ)","M(ZZ)", 200,0,1000*GeV)
pt_ZZ = ROOT.TH1F("pt(ZZ)","p_{T}(ZZ)", 200,0,1000*GeV)
numZ = ROOT.TH1F("nZ","nZ", 10,0,10)
dR = ROOT.TH1F("dR","#Delta R", 100,0,10)
ptmiss = ROOT.TH1F("ptmiss","E_{T}^{miss}",200,0,1000*GeV)

n_lep = ROOT.TH1F("nlep","n_{lep}",10,0,10)
lep_id = ROOT.TH1F("lep_id","lep_id",30,-15,-15)
lep_m = ROOT.TH1F("lep_m","m(l)",200,0,1000)
lep_pt = ROOT.TH1F("lep_pt","p_{T}(l)",200,0,1000*GeV)

#Rest of program

zz = ROOT.TChain("CollectionTree")
zz.Add("/nfs/dust/atlas/user/heim/truth5_20112017/user.sheim.mc15_13TeV.361604.PowhegPy8EG_CT10nloME_AZNLOCTEQ6L1_ZZvvll_mll4_v1_EXT0/user.sheim.12601190.EXT0._000001.DAOD_TRUTH5.TRUTH5.root")
#~ zz.Add("/nfs/dust/atlas/user/heim/truth5_20112017/user.sheim.mc15_13TeV.361604.PowhegPy8EG_CT10nloME_AZNLOCTEQ6L1_ZZvvll_mll4_v1_EXT0/user.sheim.12601190.EXT0._000002.DAOD_TRUTH5.TRUTH5.root")
#~ zz.Add("/nfs/dust/atlas/user/heim/truth5_20112017/user.sheim.mc15_13TeV.361604.PowhegPy8EG_CT10nloME_AZNLOCTEQ6L1_ZZvvll_mll4_v1_EXT0/user.sheim.12601190.EXT0._000003.DAOD_TRUTH5.TRUTH5.root")
#~ zz.Add("/nfs/dust/atlas/user/heim/truth5_20112017/user.sheim.mc15_13TeV.361604.PowhegPy8EG_CT10nloME_AZNLOCTEQ6L1_ZZvvll_mll4_v1_EXT0/user.sheim.12601190.EXT0._000004.DAOD_TRUTH5.TRUTH5.root")
#~ zz.Add("/nfs/dust/atlas/user/heim/truth5_20112017/user.sheim.mc15_13TeV.361604.PowhegPy8EG_CT10nloME_AZNLOCTEQ6L1_ZZvvll_mll4_v1_EXT0/user.sheim.12601190.EXT0._000005.DAOD_TRUTH5.TRUTH5.root")
#~ zz.Add("/nfs/dust/atlas/user/heim/truth5_20112017/user.sheim.mc15_13TeV.361604.PowhegPy8EG_CT10nloME_AZNLOCTEQ6L1_ZZvvll_mll4_v1_EXT0/user.sheim.12601190.EXT0._000006.DAOD_TRUTH5.TRUTH5.root")
#~ zz.Add("/nfs/dust/atlas/user/heim/truth5_20112017/user.sheim.mc15_13TeV.361604.PowhegPy8EG_CT10nloME_AZNLOCTEQ6L1_ZZvvll_mll4_v1_EXT0/user.sheim.12601190.EXT0._000007.DAOD_TRUTH5.TRUTH5.root")
#~ zz.Add("/nfs/dust/atlas/user/heim/truth5_20112017/user.sheim.mc15_13TeV.361604.PowhegPy8EG_CT10nloME_AZNLOCTEQ6L1_ZZvvll_mll4_v1_EXT0/user.sheim.12601190.EXT0._000008.DAOD_TRUTH5.TRUTH5.root")
#~ zz.Add("/nfs/dust/atlas/user/heim/truth5_20112017/user.sheim.mc15_13TeV.361604.PowhegPy8EG_CT10nloME_AZNLOCTEQ6L1_ZZvvll_mll4_v1_EXT0/user.sheim.12601190.EXT0._000009.DAOD_TRUTH5.TRUTH5.root")
#~ zz.Add("/nfs/dust/atlas/user/heim/truth5_20112017/user.sheim.mc15_13TeV.361604.PowhegPy8EG_CT10nloME_AZNLOCTEQ6L1_ZZvvll_mll4_v1_EXT0/user.sheim.12601190.EXT0._000010.DAOD_TRUTH5.TRUTH5.root")

t = ROOT.xAOD.MakeTransientTree(zz,ROOT.xAOD.TEvent.kBranchAccess)
total = t.GetEntries()
#~ total=10000
pas_evts=0.0
Xsec = 923.18

info = open("particle_info.txt","w")

#Selection switchboard
n_lep_switch = 0 #selects for di-lepton events
flavor_switch = 1 #selects for same flavor opp sign leptons
lep_pt_switch = 1 #selects lead and sublead lepton pt (30,20)
Z_window_switch = 1 #selects for the dilepton mass in a +/-15GeV window around Z mass
MET_switch = 0 #selects for event MET (vv pt) > 90 GeV
d_phi_switch = 0 #selects for dilepton MET phi separation > 2.7
d_R_switch = 0#selects for lepton-lepton separation R > 1.8

#Switchboard
switch = ["Lepton num switch: ","Flavor switch: ","Lepton pt switch: ","Z window switch: ","MET switch: ","dPhi switch: ","dR(ll) switch: "]

if n_lep_switch == 1: switch[0] = switch[0]+"On"
else : switch[0] = switch[0]+"Off"

if flavor_switch == 1: switch[1] = switch[1]+"On"
else : switch[1] = switch[1]+"Off"

if lep_pt_switch == 1: switch[2] = switch[2]+"On"
else : switch[2] = switch[2]+"Off"

if Z_window_switch == 1: switch[3] = switch[3]+"On"
else : switch[3] = switch[3]+"Off"

if MET_switch == 1: switch[4] = switch[4]+"On"
else : switch[4] = switch[4]+"Off"

if d_phi_switch == 1: switch[5] = switch[5]+"On"
else : switch[5] = switch[5]+"Off"

if d_R_switch == 1: switch[6] = switch[6]+"On"
else : switch[6] = switch[6]+"Off"

info.write("Switchboard\n")
for i in xrange(len(switch)):
	info.write(switch[i]+"\n")

#CutFlow counters
switch_count=[0,0,0,0,0,0,0]

MET_pass = 0
MET_fail = 0

bad_lepton = 0

#Loop over events in sample
for entry in xrange(0,total):
	
	t.GetEntry( entry )
	nZ=0
	
	lep_list=[]
	ele_list=[]
	mu_list=[]
	nu_list =[]
	Z_list = []
	jet_list=[]
	veto = 0
	
	print
	print "Event: ", t.EventInfo.eventNumber()
	#Selection of lepton candidates
	for i in xrange(t.BornLeptons.size()):
		islept = 0
		lep = t.BornLeptons[i]
		v = four_vector(lep)
		#~ if lep.status() != 1: continue
		#Candidate Selection
		if v.Pt()>7*GeV and lep.status()==23:
			if abs(lep.pdgId())==11 and abs(v.Eta())<2.47: 
				ele_list.append(lep)
				islept = 1
			#~ elif abs(lep.pdgId())==13 and abs(v.Eta())<2.5: 	
				#~ mu_list.append(lep)
				#~ islept =1 
			
			if islept == 1:
				lep_list.append(lep)
				lep_id.Fill(lep.pdgId())
		
		if abs(lep.pdgId()) in [12,14,16]:
			nu_list.append(lep)
	
	#Lepton number selection
	if len(lep_list) < 2: 
		bad_lepton+=1
		continue
	
	n_lep.Fill(len(lep_list))
	if n_lep_switch == 1:
		if len(lep_list) != 2: veto = 1  #Vetos if event does not have exactly two leptons
	
	switch_count[0]+=1
	
	sorter(ele_list)
	sorter(mu_list)
	sorter(lep_list)
	sorter(jet_list)
	
	#Same flavor opp sign selection
	tr1 = lep_list[0]
	tr2 = lep_list[1]
	if tr1.pdgId() == -(tr2.pdgId()): 
	#~ if (tr1.pdgId() == -(tr2.pdgId())) and (abs(tr1.pdgId())==11): 
		v1 = four_vector(tr1)
		v2 = four_vector(tr2)
		lep_m.Fill(v1.M())
		lep_m.Fill(v2.M())
		lep_pt.Fill(v1.Pt())
		lep_pt.Fill(v2.Pt())
		
	
	else: veto = 1
	switch_count[1]+=1
	
	#Lepton pt selection
	if lep_pt_switch == 1:
		v1 = four_vector(lep_list[0])
		v2 = four_vector(lep_list[1])
		if v1.Pt() > 30*GeV and v2.Pt() >20*GeV: pass
		else: veto = 1 
	
	switch_count[2]+=1
	
	#Z mass window selection
	if Z_window_switch == 1:
		v1 = four_vector(lep_list[0])
		v2 = four_vector(lep_list[1])
		lep2 = v1+v2
		M = lep2.M()

		if M>76*GeV and M<106*GeV: 
			ll_mass.Fill(M)
			ll_pt.Fill(lep2.Pt())
			ll_eta.Fill(lep2.Rapidity())
			
		else : veto = 1
	
	switch_count[3]+=1
	
	netMET = ROOT.TLorentzVector()
		
	sorter(nu_list)
	
	ptmiss.Fill(netMET.Pt())
	MET=ROOT.TLorentzVector()
	
	#MET selection
	if len(nu_list)<2: veto = 1
	tr1 = nu_list[0]
	tr2 = nu_list[1]
	if tr1.pdgId() == -(tr2.pdgId()):
		v1=four_vector(tr1)
		v2=four_vector(tr2)
		MET = v1+v2
		vv_mass.Fill(MET.M())
		vv_pt.Fill(MET.Pt())
		vv_eta.Fill(MET.Rapidity())
		MET_pass+=1
	else: veto = 1
	
	if veto == 0 :
		pas_evts+=1
	else:
	#Dumps particle info
		print >>info,"\nEvent: "+str(t.EventInfo.eventNumber())+"\t Run: "+str(t.EventInfo.runNumber())
		for i in xrange(t.BornLeptons.size()):
			tr = t.BornLeptons[i]
			print >>info, "%d \t %+d \t %d \t %d \t %0.4f \t %0.4f \t %0.4f \t %0.4f" %(i, tr.pdgId(),tr.status(),tr.nParents(),tr.pt()/GeV,tr.eta(),tr.phi(),tr.m()/GeV)
			
		for i in xrange(len(nu_list)):
			tr = nu_list[i]
			print >>info, "%d \t %+d \t %d \t %d \t %0.4f \t %0.4f \t %0.4f \t %0.4f" %(i, tr.pdgId(),tr.status(),tr.nParents(),tr.pt()/GeV,tr.eta(),tr.phi(),tr.m()/GeV)
	
# Back to main body
file1.Write()
print "Event fraction passing cuts: ",pas_evts,"/",total
print "Effective cross section: ", pas_evts/total*Xsec, " fb"
print "Bad events with < 2 leptons: ",bad_lepton
ROOT.xAOD.ClearTransientTrees()
print("Done!!")
info.close()
