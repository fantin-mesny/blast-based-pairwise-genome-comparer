#!/usr/bin/python2
#developed in august 2017, by Fantin Mesny
#contact: fantin.mesny AT icloud.com


import argparse
import sys
import subprocess
import os
import glob
import csv
import matplotlib.pyplot as plt
from Bio import SeqIO
import networkx as nx
import datetime



def get_params(argv, reqf, reqo):
	parser = argparse.ArgumentParser(description='A program to compare two similar genomes. Specify either two PROKKA-output folders (-folders) OR pairs of transcripts (-rna), proteins (-prot), genomes (-genomes) and gff (-gff) files.')
	parser.add_argument('-folders', '--folders', help='Two comma-separated folders including transcripts (.ffn), proteins (.faa) and genome (.fna) fasta files', required=reqf) 
	parser.add_argument('-genomes', '--genomes', help='Two comma-separated genomes (fasta-files)', required=reqo)
	parser.add_argument('-prot', '--prot', help='Two comma-separated predicted proteins multi-fasta files', required=reqo)
	parser.add_argument('-rna', '--rna', help='Two comma-separated predicted transcripts multi-fasta files', required=reqo)
	parser.add_argument('-gff', '--gff', help='Two comma-separated .gff annotation files', required=reqo)
	parser.add_argument('-o', '--o', help='Output folder', required=False)
	parser.add_argument('-minLrap','--minLrap', help='minLrap significance threshold for blast hits', default=0.7)
	parser.add_argument('-pident','--pident', help='Identity percentage significance threshold for blast hits', default=70.)
	parser.add_argument('-graph','--graph', help='Ouput type for graphical representations. "img" (default) or "show" to display using Python.', default='img')
	parser.add_argument('-cog','--cog', help="Perform COG rpsblasts ? 'yes'(default) or 'no'", default='no')
	parser.add_argument('-pfam','--pfam', help="Perform Pfam rpsblasts and get protein-GO associations ? 'yes'(default) or 'no'", default='no')
	parser.add_argument('-paralog_pid','--paralog_pid', help='Identity percentage significance threshold to consider proteins as homologous (default: 60)', default=60)
	a = parser.parse_args()
	return a

def terminal(cmd, Dir):
	#Dir=a.o
	p = subprocess.Popen("cd " + Dir + " ; " + cmd, shell=True, stdout = subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
	p.wait()

def getFileWithExtension(directory, extension):
	os.chdir(directory)
	i=0
	for fil in glob.glob("*." + extension):
		filename = directory + fil
		i+=1
	if i==1:
		return filename
	elif i==0:
		print("ERROR: no ." + extension + " file found in directory " + directory)
	elif i>1:
		print("ERROR: multiple ." + extension + " files in given directory " + directory)

def getHitsTable(q1,q2,b1,b2,o):
	protsq1=SeqIO.to_dict(SeqIO.parse(q1, "fasta"))
	protsq2=SeqIO.to_dict(SeqIO.parse(q2, "fasta"))
	corresp21=[line for line in csv.DictReader(open(b1,"r"), delimiter='\t', fieldnames=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])]
	corresp12=[line for line in csv.DictReader(open(b2,"r"), delimiter='\t', fieldnames=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])]
	with open(o, 'w+') as csvfile:
   		fieldnames = ['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore','minLrap','maxLrap','reciprocity']
   		output = csv.DictWriter(csvfile, delimiter='\t', fieldnames=fieldnames)
		for hit in corresp21:
			align=float(hit['length'])
			if hit['qseqid'] in protsq1:
				Len1=float(len(protsq1[hit['qseqid']].seq))
				Len2=float(len(protsq2[hit['sseqid']].seq))
			else:
				Len1=float(len(protsq1[hit['sseqid']].seq))
				Len2=float(len(protsq2[hit['qseqid']].seq))			
			if Len1>Len2:
				hit['minLrap']=round(align/Len2,3)
				hit['maxLrap']=round(align/Len1,3)
			else:
				hit['minLrap']=round(align/Len1,3)
				hit['maxLrap']=round(align/Len2,3)
			if any([hit12['qseqid']==hit['sseqid'] and hit12['sseqid']==hit['qseqid'] for hit12 in corresp12]):
				hit['reciprocity']='Reciprocal'
			else:
				hit['reciprocity']='.'
			output.writerow({x:hit[x] for x in fieldnames})

def keepCDS(p1,p2,tx1,tx2):
	pr1=[x.id for x in SeqIO.parse(p1, "fasta")]
	pr2=[x.id for x in SeqIO.parse(p2, "fasta")]
	corresp21=[line for line in csv.DictReader(open(tx1,"r"), delimiter='\t', fieldnames=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])]
	corresp12=[line for line in csv.DictReader(open(tx2,"r"), delimiter='\t', fieldnames=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])]
	with open(tx1, 'w') as outp1:
   		fieldnames = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
   		output = csv.DictWriter(outp1, delimiter='\t', fieldnames=fieldnames)
		removed=0
		for hit in corresp21:
			if hit['qseqid'] in pr2 and hit['sseqid'] in pr1:
				output.writerow({x:hit[x] for x in fieldnames})
			else:
				removed+=1
	print("_ "+str(removed)+" hits have been removed from tblastx 1")
	with open(tx2, 'w') as outp2:
   		fieldnames = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
   		output = csv.DictWriter(outp2, delimiter='\t', fieldnames=fieldnames)
		removed=0
		for hit in corresp12:
			if hit['qseqid'] in pr1 and hit['sseqid'] in pr2:
				output.writerow({x:hit[x] for x in fieldnames})
			else:
				removed+=1
	print("_ "+str(removed)+" hits have been removed from tblastx 2")


def CDStable(hitstableblastp,hitstabletblastx,gff,minLrap,pident,o):
	hitTableFields=['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore', 'minLrap', 'maxLrap', 'reciprocity']
	gffFields=['strain', 'program','type', 'cds_start','cds_end', 'dot', 'sense', 'o', 'annot']
	hitslistP=[line for line in csv.DictReader(open(hitstableblastp,"r"), delimiter='\t', fieldnames=hitTableFields)]
	hitslistT=[line for line in csv.DictReader(open(hitstabletblastx,"r"), delimiter='\t', fieldnames=hitTableFields)]
	cdslist=[line for line in csv.DictReader(open(gff,"r"), delimiter='\t', fieldnames=gffFields) if '#' not in line['strain'] and line['type']]
	with open(o,'w+') as output:
		nbcds=0
		for cds in cdslist: #chaque cds prise une par une
			a=0
			cdsid=cds['annot'].split(';')[0][3:]
			if cds['type']=="CDS":
				nbcds+=1
				for hit in hitslistP: #cds recherchee dans blastp
					if hit['qseqid']==cdsid: #and float(hit['minLrap'])>=minLrap and float(hit['pident'])>pident:
						for hitT in hitslistT: #cds trouvee dans le blastp, maintenant recherchee dans tblastx
							if hitT['qseqid']==cdsid and hitT['sseqid']==hit['sseqid']: 
								output.write(cdsid+'\t'+cds['cds_start']+'\t'+cds['cds_end']+'\t'+hit['sseqid']+'\t'+hit['pident']+"\t"+hit['minLrap']+'\t'+hit['reciprocity']+'\t'+hitT['pident']+"\t"+hitT['minLrap']+'\t'+hitT['reciprocity']+'\n')
								a=1
						if a==0:#cds pas trouvee dans le tblastx mais trouve dans blastp
							output.write(cdsid+'\t'+cds['cds_start']+'\t'+cds['cds_end']+'\t'+hit['sseqid']+'\t'+hit['pident']+"\t"+hit['minLrap']+'\t'+hit['reciprocity']+'\t0.00\t0.00\t.\n')
							a=1
				if a==0: #cds pas trouvee dans le blastp, recherche dans le tblastx
					for hitT in hitslistT:
							if hitT['qseqid']==cdsid: #and float(hitT['minLrap'])>minLrap and float(hitT['pident'])>pident:
								output.write(cdsid+'\t'+cds['cds_start']+'\t'+cds['cds_end']+'\t'+hitT['sseqid']+'\t0.00\t0.00\t.\t'+hitT['pident']+'\t'+hitT['minLrap']+'\t'+hitT['reciprocity']+'\n')
								a=1
				if a==0: #cds sans hit dans le tblastx ni dans le blastp
					output.write(cdsid+'\t'+cds['cds_start']+'\t'+cds['cds_end']+'\t.\t0.00\t0.00\t.\t0.00\t0.00\t.\n')


def reciprocalBestMatches(cdslist,minLrap,pident,o):
	cdslist=[line for line in csv.DictReader(open(cdslist,"r"), delimiter='\t', fieldnames=['qseqid','CDSstart','CDSend','sseqid','ppident','pmin', 'preciprocity', 'tpident','tmin', 'treciprocity'])]
	subjects=[]
	with open(o,'w+') as output:
		for cds in cdslist:
			if 'eciprocal' in cds['preciprocity'] and float(cds['ppident'])>=pident and float(cds['pmin'])>=minLrap:
				subjects.append(cds['sseqid'])
				output.write(cds['qseqid']+"\t"+cds['sseqid']+"\n")
		for cds in cdslist:
			if 'eciprocal' in cds['treciprocity'] and float(cds['tpident'])>=pident and float(cds['tmin'])>=minLrap and cds['sseqid'] not in subjects:
				output.write(cds['qseqid']+"\t"+cds['sseqid']+"\n")


def CDSgraph(name1,ht1,name2,ht2,o,graph):
	i1=[line for line in csv.DictReader(open(ht1,"r"), delimiter='\t', fieldnames=['cds_id', 'cds_start', 'cds_stop', 'corresp', 'pident', 'minLrap', 'reciprocity'])]
	fig1=plt.figure()
	ax = fig1.add_subplot(2,1,1)
	ax.title.set_text(name1)
	x1=[a['cds_start'] for a in i1 if 'ecipro' in a['reciprocity']]
	x2=[a['cds_start'] for a in i1 if 'ecipro' not in a['reciprocity']]
	ax.plot(x1,[a['pident'] for a in i1 if 'ecipro' in a['reciprocity']], '.', color='red', alpha=0.5)
	ax.plot(x2,[a['pident'] for a in i1 if 'ecipro' not in a['reciprocity']], '.', color='blue', alpha=0.5)
	plt.xlabel("CDS")
	plt.ylabel("%Identity")
	ax2 = fig1.add_subplot(2,1,2)
	i2=[line for line in csv.DictReader(open(ht2,"r"), delimiter='\t', fieldnames=['cds_id', 'cds_start', 'cds_stop', 'corresp', 'pident', 'minLrap', 'reciprocity'])]
	x12=[a['cds_start'] for a in i2 if 'ecipro' in a['reciprocity']]
	x22=[a['cds_start'] for a in i2 if 'ecipro' not in a['reciprocity']]
	ax2.plot(x12,[a['pident'] for a in i2 if 'ecipro' in a['reciprocity']], '.', color='red', alpha=0.5)
	ax2.plot(x22,[a['pident'] for a in i2 if 'ecipro' not in a['reciprocity']], '.', color='blue', alpha=0.5)
	ax2.title.set_text(name2)
	fig1.subplots_adjust(hspace=.5)
	plt.xlabel("CDS")
	plt.ylabel("%Identity")
	if "show" in graph:
		plt.show()
	else:
		plt.savefig(o, dpi=500)

def getSignificantCogs(rpsblast,thresholds,out,outname):
	thresh={line.split('\t')[0]:line.split('\t')[2] for line in open(thresholds,'r').readlines() if line.split('\t')[1][:3]=="COG"}
	cddtocog={line.split('\t')[0]:line.split('\t')[1] for line in open(thresholds,'r').readlines() if line.split('\t')[1][:3]=="COG"}
	rps=[line for line in csv.DictReader(open(rpsblast,"r"), delimiter='\t', fieldnames=['qseqid', 'sseqid', 'qstart', 'qend', 'sstart', 'send', 'qlen', 'slen', 'bitscore', 'pident', 'nident', 'mismatch', 'evalue'])]
	dicog={}
	with open(out+outname+'.cogperprot.csv','w+') as output:
		for a in rps:
			cddid=a['sseqid'].split('|')[2]
			if float(thresh[cddid])<=float(a['bitscore']):
				output.write(a['qseqid']+'\t'+cddtocog[cddid]+'\t'+cddid+'\t'+a['bitscore']+'\t'+a['evalue']+'\n')
				if cddid not in dicog:
					dicog[cddid]=[]
					dicog[cddid].append(a['qseqid'])
				elif a['qseqid'] not in dicog[cddid]:
					dicog[cddid].append(a['qseqid'])

def getSignificantPfam(rpsblast,thresholds,pfam2go,out,outname):
	thresh={line.split('\t')[0]:line.split('\t')[2] for line in open(thresholds,'r').readlines() if line.split('\t')[1][:4]=="pfam"}
	cddtopfam={line.split('\t')[0]:line.split('\t')[1] for line in open(thresholds,'r').readlines() if line.split('\t')[1][:4]=="pfam"}
	pfamtogo={}
	for line in open(pfam2go, "r").readlines():
		if line[0] !='!':
			if line.split(' ')[0][7:] not in pfamtogo:
				pfamtogo[line.split(' ')[0][7:]]=[]
			pfamtogo[line.split(' ')[0][7:]].append(line.split(' ')[-1][:-1])
	rps=[line for line in csv.DictReader(open(rpsblast,"r"), delimiter='\t', fieldnames=['qseqid', 'sseqid', 'qstart', 'qend', 'sstart', 'send', 'qlen', 'slen', 'bitscore', 'pident', 'nident', 'mismatch', 'evalue'])]
	dicog={}
	with open(out+outname+'.pfamperprot.csv','w+') as output:
		for a in rps:
			cddid=a['sseqid'].split('|')[2]
			if float(thresh[cddid])<=float(a['bitscore']):
				output.write(a['qseqid']+'\t'+cddtopfam[cddid]+'\t'+cddid+'\t'+a['bitscore']+'\t'+a['evalue']+'\n')
				if cddid not in dicog:
					dicog[cddid]=[]
					dicog[cddid].append(a['qseqid'])
				elif a['qseqid'] not in dicog[cddid]:
					dicog[cddid].append(a['qseqid'])
	with open(out+outname+'.goperprot.csv','w+') as output:
		for a in rps:
			cddid=a['sseqid'].split('|')[2]
			if float(thresh[cddid])<=float(a['bitscore']):
				if cddtopfam[cddid][4:] in pfamtogo:
					output.write(a['qseqid']+'\t'+';'.join(pfamtogo[cddtopfam[cddid][4:]])+'\n')
				else:
					pass
				if cddid not in dicog:
					dicog[cddid]=[]
					dicog[cddid].append(a['qseqid'])
				elif a['qseqid'] not in dicog[cddid]:
					dicog[cddid].append(a['qseqid'])

def paralogFinder(makeblastdb,a,p1,p2,out):
	pid=60
	id1=[seq for seq in SeqIO.to_dict(SeqIO.parse(p1, "fasta"))]
	id2=[seq for seq in SeqIO.to_dict(SeqIO.parse(p2, "fasta"))]
	if makeblastdb=="yes":
		terminal("makeblastdb -in "+ p1 +" -dbtype 'prot' -title 'prot1_db' -out 'prot1_db'", out)
		terminal("makeblastdb -in "+ p2 +" -dbtype 'prot' -title 'prot2_db' -out 'prot2_db'", out)
	terminal("blastp -db 'prot1_db' -query '" + p2 + "' -out '2on1.blastp_forpf' -outfmt '6' -evalue 0.001 -max_hsps 1", out)
	terminal("blastp -db 'prot2_db' -query '" + p1 + "' -out '1on2.blastp_forpf' -outfmt '6' -evalue 0.001 -max_hsps 1", out)
	terminal("blastp -db 'prot2_db' -query '" + p2 + "' -out '2on2.blastp_forpf' -outfmt '6' -evalue 0.001 -max_hsps 1", out)
	terminal("blastp -db 'prot1_db' -query '" + p1 + "' -out '1on1.blastp_forpf' -outfmt '6' -evalue 0.001 -max_hsps 1", out)
	print("_Running blasts...")
	blast21=[line for line in csv.DictReader(open(out+"2on1.blastp_forpf","r"), delimiter='\t', fieldnames=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])]
	blast12=[line for line in csv.DictReader(open(out+"1on2.blastp_forpf","r"), delimiter='\t', fieldnames=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])]
	blast11=[line for line in csv.DictReader(open(out+"1on1.blastp_forpf","r"), delimiter='\t', fieldnames=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])]
	blast22=[line for line in csv.DictReader(open(out+"2on2.blastp_forpf","r"), delimiter='\t', fieldnames=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])]
	out=out+'paralogs/'
	queries21={}
	queries12={}
	queries11={}
	queries22={}
	print("_Analysing blasts results...")
	for line in blast21:
		if line['qseqid'] not in queries21:
			queries21[line['qseqid']]=[]
		if float(line['pident'])>pid: 
			queries21[line['qseqid']].append((line['sseqid'],line['pident']))
	queries21={a:queries21[a] for a in queries21 if len(queries21[a])>1}
	for line in blast12:
		if line['qseqid'] not in queries12:
			queries12[line['qseqid']]=[]
		if float(line['pident'])>pid: 
			queries12[line['qseqid']].append((line['sseqid'],line['pident']))
	queries12={a:queries12[a] for a in queries12 if len(queries12[a])>1}
	for line in blast11:
		if line['qseqid'] not in queries11:
			queries11[line['qseqid']]=[]
		if float(line['pident'])>pid and line['sseqid']!=line['qseqid']: 
			queries11[line['qseqid']].append((line['sseqid'],line['pident']))
	queries11={a:queries11[a] for a in queries11 if len(queries11[a])>1}
	for line in blast22:
		if line['qseqid'] not in queries22:
			queries22[line['qseqid']]=[]
		if float(line['pident'])>pid and line['sseqid']!=line['qseqid']: 
			queries22[line['qseqid']].append((line['sseqid'],line['pident']))
	queries22={a:queries22[a] for a in queries22 if len(queries22[a])>1}
	print("_Assembling paralog families...")
	listedegraphes=[]
	listedesnoeuds=[]
	with open(out+'paralogfamilies.csv', 'w+') as outcsv:
		output = csv.DictWriter(outcsv, delimiter='\t', fieldnames=['family_number','ids_from_1','ids_from_2'])
		i=0
		for prot in queries21:
			l=[]
			idsfrom1=[]
			idsfrom2=[]
			idsfrom2.append(prot)
			for s in queries21[prot]:
				l.append((prot,s[0],s[1]))
				idsfrom1.append(s[0])
				if s[0] in queries11:
					for s4 in queries11[s[0]]:
						l.append((s[0],s4[0],s4[1]))
						if s4[0] not in idsfrom1:
							idsfrom1.append(s4[0])
				if s[0] in queries12:
					for s2 in queries12[s[0]]:
						l.append((s[0],s2[0],s2[1]))
						if s2[0] not in idsfrom2:
							idsfrom2.append(s2[0])
						if s2[0] in queries22:
							for s3 in queries22[s2[0]]:
								l.append((s2[0],s3[0],s3[1]))
								if s2[0] not in idsfrom2:
									idsfrom2.append(s3[0])
			if not any([set(idsfrom1+idsfrom2).issubset(nodelist) for nodelist in listedesnoeuds]):
				i+=1
				listedegraphes.append(l)
				listedesnoeuds.append(idsfrom1+idsfrom2)
				output.writerow({'family_number':i,'ids_from_1':','.join(idsfrom1),'ids_from_2':','.join(idsfrom2)})
	print("_Building networks...")
	i=0
	for graph in listedegraphes:
		i+=1
		formatted_list = [(node[0], node[1], {"weight":float(node[2])/100., "label":node[2]}) for node in graph]
		plt.figure()
		G = nx.Graph()
		G.add_edges_from(formatted_list)
		edges,weights = zip(*nx.get_edge_attributes(G,'weight').items())
		color_map=[]
		for node in G:
			if node in id1:
				color_map.append('blue')
			else: 
				color_map.append('green')
		nx.draw(G, node_color=color_map, edgelist=edges, arrows=True, edge_color=weights, edge_cmap=plt.cm.Blues, linewidths=0, alpha=0.5)
		plt.savefig(out+str(i)+'.png')
		nx.write_graphml(G, out+str(i)+".graphml")

def corresp2go(goperprot1,goperprot2,out,corresp):
	dic1={}
	dic2={}
	with open(goperprot1, "r") as go1:
		for line in go1.readlines():
			dic1[line.split('\t')[0]]=line[:-1].split('\t')[1]
	with open(goperprot2, "r") as go2:
		for line in go2.readlines():
			dic2[line.split('\t')[0]]=line[:-1].split('\t')[1]
	with open(corresp, "r") as corresp:
		with open(out+"corresp2go.csv","w+") as c2go:
			for line in corresp.readlines():
				if line.split('\t')[1][:-1] in dic1 and line.split('\t')[0] in dic2:
					c2go.write(line.split('\t')[1][:-1]+'-'+line.split('\t')[0]+'\t'+';'.join(list(set(dic1[line.split('\t')[1][:-1]].split(';')+dic2[line.split('\t')[0]].split(';'))))+'\n')
				else:
					c2go.write(line.split('\t')[1][:-1]+'-'+line.split('\t')[0]+'\n')
			
	


if __name__ == '__main__':
	if "-folders" in sys.argv[1:]:
		reqf=True
		reqo=False
	elif "-genomes" in sys.argv[1:]:
		reqf=False
		reqo=True
	elif "-h" in sys.argv[1:]:
		reqf=False
		reqo=False
	else:
		reqf=False
		reqo=False
		print ("ERROR: Please specify either PROKKA folders (-folders) or pairs of transcripts (-rna), proteins (-prot), genomes (-genomes) and gff (-gff) files.")
	a = get_params(sys.argv[1:],reqf,reqo)
	appdir= str(os.path.dirname(os.path.realpath(__file__)))+"/"
	cwd = os.getcwd()
	now = datetime.datetime.now()
	date=str(now.year)+str(now.month)+str(now.day)+"_"+str(now.hour)+str(now.minute)+str(now.second)
	if "~" in str(a.folders) or "~" in str(a.genomes) or "~" in str(a.gff) or "~" in str(a.rna) or "~" in str(a.prot):
		print("ERROR: '~' character is not allowed")
	if "../" in str(a.folders) or "../" in str(a.genomes) or "../" in str(a.gff) or "../" in str(a.rna) or "../" in str(a.prot):
		print("ERROR: '../' directories are not allowed, please use the absolute directory")
	if a.o:
		out=a.o
		if out[0]!='/':
			out=cwd+"/"+out
		if out[-1]!='/':
			out=out+"/"
	else:
		terminal("mkdir "+date, appdir)
		out=appdir+date+"/"
	
	if a.folders:
		folders=a.folders
		path1=folders.split(",")[0]
		path2=folders.split(",")[1]
		if path1[0]!='/':
			out=cwd+"/"+out
		if path1[-1] != '/':
			path1=path1+'/'
		if path2[0]!='/':
			out=cwd+"/"+out
		if path2[-1] != '/':
			path2=path2+'/'
		genome1=getFileWithExtension(path1, "fna")
		genome2=getFileWithExtension(path2, "fna")
		prot1=getFileWithExtension(path1, "faa")
		prot2=getFileWithExtension(path2, "faa")
		rna1=getFileWithExtension(path1, "ffn")
		rna2=getFileWithExtension(path2, "ffn")
		gff1=getFileWithExtension(path1, "gff")
		gff2=getFileWithExtension(path2, "gff")
	else:
		genome1=a.genomes.split(",")[0]
		genome2=a.genomes.split(",")[1]
		prot1=a.prot.split(",")[0]
		prot2=a.prot.split(",")[1]
		rna1=a.rna.split(",")[0]
		rna2=a.rna.split(",")[1]
		gff1=a.gff.split(",")[0]
		gff2=a.gff.split(",")[1]
	minLrap=float(a.minLrap)
	pident=float(a.pident)
	name1=open(genome1, "r").readlines()[0][1:-1]
	name2=open(genome2, "r").readlines()[0][1:-1]
	cog=str(a.cog)
	pfam=str(a.pfam)
	graph=str(a.graph)
	if graph!="img" or graph!="show":
		print("ERROR: incorrect argument for -graph, images will be made")
	print(prot1,prot2)

	print('\nArguments correctly parsed')
	print("Output folder: "+out)
	print("Creating blastp database on first protein file...")
	terminal("makeblastdb -in "+ prot1 +" -dbtype 'prot' -title 'prot1_db' -out 'prot1_db'", out)
	print("Performing first blastp (protein file 2 on protein file 1)...")
	terminal("blastp -db 'prot1_db' -query '" + prot2 + "' -out '2on1.blastp' -outfmt '6' -max_target_seqs 1 -evalue 0.001 -max_hsps 1", out)
	print("Creating blastp database on second protein file...")
	terminal("makeblastdb -in "+ prot2 +" -dbtype 'prot' -title 'prot2_db' -out 'prot2_db'", out)
	print("Performing second blastp (protein file 1 on protein file 2)...")
	terminal("blastp -db 'prot2_db' -query '" + prot1 + "' -out '1on2.blastp' -outfmt '6' -max_target_seqs 1 -evalue 0.001 -max_hsps 1", out)
	print("Running paralogfinder...")
	terminal("mkdir paralogs", out)
	paralogFinder("no",a,prot1,prot2,out)
	terminal("rm *_forpf", out)
	print("Removing blast-db files")
	terminal("rm prot1_db* ; rm prot2_db*", out)
	print("Building 'blastp2on1_hitsTable.csv' file...")
	getHitsTable(prot1,prot2,out+"2on1.blastp",out+"1on2.blastp",out+'blastp2on1_hitsTable.csv')
	print("Building 'blastp1on2_hitsTable.csv' file...")
	getHitsTable(prot2,prot1,out+"1on2.blastp",out+"2on1.blastp",out+'blastp1on2_hitsTable.csv')

	print("Creating blast database on first transcript file...")
	terminal("makeblastdb -in "+ rna1 +" -dbtype 'nucl' -title 'rna1_db' -out 'rna1_db'", out)
	print("Performing first tblastx (transcript file 2 on transcript file 1)...")
	terminal("tblastx -db 'rna1_db' -query " + rna2 + " -out '2on1.tblastx' -outfmt '6' -max_target_seqs 1 -seg 'no' -max_hsps 1", out)
	print("Creating blast database on second protein file...")
	terminal("makeblastdb -in "+ rna2 +" -dbtype 'nucl' -title 'rna2_db' -out 'rna2_db'", out)
	print("Performing second tblastx (transcript file 1 on transcript file 2)...")
	terminal("tblastx -db 'rna2_db' -query '" + rna1 + "' -out '1on2.tblastx' -outfmt '6' -max_target_seqs 1 -seg 'no' -max_hsps 1", out)
	print("Removing blast-db files")
	terminal("rm rna1_db* ; rm rna2_db*", out)
	print("Removing non-coding transcripts data from tblastx files...")
	keepCDS(prot1,prot2,out+"2on1.tblastx",out+"1on2.tblastx")
	print("Building 'tblastx2on1_hitsTable.csv' file...")
	getHitsTable(prot1,prot2,out+"2on1.tblastx",out+"1on2.tblastx",out+'tblastx2on1_hitsTable.csv')
	print("Building 'tblastx1on2_hitsTable.csv' file...")
	getHitsTable(prot2,prot1,out+"1on2.tblastx",out+"2on1.tblastx",out+'tblastx1on2_hitsTable.csv')
	print("Building 'CDStable2on1.csv' file...")
	CDStable(out+"blastp2on1_hitsTable.csv",out+"tblastx2on1_hitsTable.csv",gff2,minLrap,pident,out+'CDStable2on1.csv')
	print("Building 'CDStable1on2.csv' file...")
	CDStable(out+"blastp1on2_hitsTable.csv",out+"tblastx1on2_hitsTable.csv",gff1,minLrap,pident,out+'CDStable1on2.csv')
	print("Extracting a list of reciprocal best matches...")
	reciprocalBestMatches(out+'CDStable2on1.csv',minLrap,pident,out+'ReciprocalBestMatches.csv')
	print("Building figures...")
	CDSgraph(name1,out+'CDStable1on2.csv',name2,out+'CDStable2on1.csv',out+'Genomes_similarity_graph.png',graph)

## Uncomment the following lines to get functional data: gene-GOs, gene-Pfam and gene-COG associations
## You need to create a "COG" folder including a COG blast database and a "Pfam" one containing a Pfam blast database into this script directory
## You also need to put a csv file with bitscore threshold (on column 2) for each pfam/cog number. Modify the getSignificant* functions to make them correctly parse this file.
 
#	if "no" not in cog:
#		print("Performing RPSblasts on COG database...")
#		terminal("rpsblast -query "+prot1+" -db "+appdir+"/COG/Cog -evalue 0.001 -outfmt '6 qseqid sseqid qstart qend sstart send qlen slen bitscore pident nident mismatch evalue' -out "+out+"1.cog.rpsblast",out)
#		getSignificantCogs(out+"1.cog.rpsblast",appdir+"thresholds.txt",out,"1")
#		terminal("rpsblast -query "+prot2+" -db "+appdir+"/COG/Cog -evalue 0.001 -outfmt '6 qseqid sseqid qstart qend sstart send qlen slen bitscore pident nident mismatch evalue' -out "+out+"2.cog.rpsblast",out)
#		getSignificantCogs(out+"2.cog.rpsblast",appdir+"thresholds.txt",out,"2")
#	if "no" not in pfam:
#		print("Performing RPSblasts on Pfam database...")
#		terminal("rpsblast -query "+prot1+" -db "+appdir+"Pfam/Pfam -evalue 0.001 -outfmt '6 qseqid sseqid qstart qend sstart send qlen slen bitscore pident nident mismatch evalue' -out "+out+"1.pfam.rpsblast",out)
#		getSignificantPfam(out+"1.pfam.rpsblast",appdir+"thresholds.txt",appdir+"Pfam/pfam2go",out,"1")
#		terminal("rpsblast -query "+prot2+" -db "+appdir+"Pfam/Pfam -evalue 0.001 -outfmt '6 qseqid sseqid qstart qend sstart send qlen slen bitscore pident nident mismatch evalue' -out "+out+"2.pfam.rpsblast",out)
#		getSignificantPfam(out+"2.pfam.rpsblast",appdir+"thresholds.txt",appdir+"Pfam/pfam2go",out,"2")
#		#ecriture du fichier corresp2go necessaire aux analyses statistiques GO avec topGO
#		corresp2go(out+'1.goperprot.csv',out+'2.goperprot.csv',out,out+'ReciprocalBestMatches.csv')
	
