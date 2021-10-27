#!/usr/bin/env python3

import os, sys, subprocess, argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


class Pesti():
		def __init__(self,path,outdir,threads,refs):
			self.path = path
			self.outdir = outdir
			self.threads = threads
			self.refs = refs
		
		def run_fastqc(self,path,outdir):
			print("\n")
			if not os.path.exists(os.path.join(outdir,"fastqc")):
				os.mkdir(os.path.join(outdir,"fastqc"))
			for file in os.listdir(path):
				(base,ext)=os.path.splitext(file)
				if ext == ".fastq":
					print(f"Running fastqc on {os.path.join(path,file)}")
					subprocess.run(['fastqc','-o',os.path.join(outdir,"fastqc"), os.path.join(path,file)],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

		def trimming(self, path, outdir):
			adapters="../BBR_subsample_virus/adapters.fa"
			print("\n")
			if not os.path.exists(os.path.join(outdir,"1_trimming")):
				os.mkdir(os.path.join(outdir,"1_trimming"))
			for file in os.listdir(path):
				(base,ext)=os.path.splitext(file)
				if ext == ".fastq":
					print(f"Running bbduk on {os.path.join(path,file)}")
					subprocess.run(['bbduk.sh','-Xmx2g',f'in1={os.path.join(path,file)}',f'out1={os.path.join(outdir,"1_trimming/",base+".trim.fastq")}',f'stats={os.path.join(outdir,"1_trimming",base+".stats.txt")}',f'ref={adapters}','ktrim=r','k=23','mink=11','hdist=1','qtrim=rl','trimq=10','maq=10','minlen=50','threads=48','tpe','tbo'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

		def map2ref(self, path, outdir, threads,refs):
			print("\n")
			if not os.path.exists(os.path.join(outdir,"2_map2ref")):
				os.mkdir(os.path.join(outdir,"2_map2ref"))
			subprocess.run(['bwa','index','-p','PIPi_index_temp',refs],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
			for file in os.listdir(os.path.join(outdir,"1_trimming")):
				(base,ext)=os.path.splitext(file)
				if ext == ".fastq":
					print(f"Mapping reads from {file} to data/pestivirus_full_genomes.fasta")
					sam=subprocess.Popen(['bwa','mem','-t',threads,'PIPi_index_temp',os.path.join(outdir,"1_trimming",file)],stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
					sort=subprocess.Popen(['samtools','sort','-'],stdin=sam.stdout,stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
					subprocess.run(['samtools','coverage', '-o',os.path.join(outdir,"2_map2ref",base+".sort.cov"),'-'],stdin=sort.stdout)
					print(f"Calculating coverage for {file}\n")


		def output_graph(self,dataframe,outdir,base):
			if not os.path.exists(os.path.join(outdir,"graphs")):
				os.mkdir(os.path.join(outdir,"graphs"))
			subset=dataframe[["#rname","Percentage"]]
			other = 100 - subset['Percentage'].sum()
			df = pd.DataFrame([["Other", other]], columns=["#rname","Percentage"])
			subset=subset.append(df)
			subset.set_index('#rname', inplace=True)
			sizes=subset['Percentage']
			# I'm not sure about this function inside here... but it works, so I'm sticking with it.
			def absolute_value(val):
				a= np.round(val/100.*sizes.sum(), 1)
				return a
			plot=subset.plot.pie(y='Percentage', figsize=(20,20), labeldistance=0.5, autopct=absolute_value)
			plot.get_legend().remove()
			plt.savefig(os.path.join(outdir,'graphs',base+'.png'))

		def parse_cov(self,outdir):
			for file in os.listdir(os.path.join(outdir,"2_map2ref")):
				(base,ext)=os.path.splitext(file)
				if ext == ".cov":
					df=pd.read_csv(os.path.join(outdir,"2_map2ref",file), sep="\t")
					total = df['numreads'].sum()
					percentage=df["numreads"] / total *100
					df["Percentage"]=percentage
					sorted=df.loc[(df['Percentage'] >= 5)].sort_values(by="Percentage", ascending=False).round({'Percentage':2})
					sorted.to_csv(os.path.join(outdir,'2_map2ref',base+'cov_results.csv'))
					self.output_graph(sorted,outdir,base)

		def run(self):
			if not os.path.exists(self.outdir):
				os.mkdir(self.outdir)
			else:
				answer=""
				while answer.lower() != "y" or answer.lower() != "n":
					answer=input("Previous output detected, overwite? (y/n): ")
					if answer == "n":
						print("Exiting...")
						sys.exit()
					elif answer != "y":
						print("Please chose either 'y' or 'n'")
					elif answer == "y":
						break
			self.run_fastqc(self.path, self.outdir)
			self.trimming(self.path,self.outdir)
			self.map2ref(self.path, self.outdir, self.threads,self.refs)
			self.parse_cov(self.outdir)
			#TBC
			#self.IRMA 
			#self.assembly

def main():
	parser = argparse.ArgumentParser(description="Pestivirus Identification Pipeline", formatter_class=argparse.RawTextHelpFormatter, add_help=False)
	required = parser.add_argument_group('Required Arguments')
	optional = parser.add_argument_group('Optional Arguments')
	

	required.add_argument('-p', '--path', type=str, required=True, help="Folder of fastq files from Ion torrent")
	required.add_argument("-o", "--outdir",type=str, required=True, help="Output directory")
	required.add_argument("-r", "--reference", type=str, required=True, help="Location of reference Pestivirus sequences")
	optional.add_argument("-t", "--threads", type=str, required=False, default=8, help="Number of threads to use. Default 8")
	optional.add_argument("-h", "--help", action="help", help="show this help message and exit")


	args=parser.parse_args()
	path=os.path.normpath(args.path)
	outdir=args.outdir
	refs=args.reference
	threads=args.threads

	job=Pesti(path,outdir,threads,refs)
	job.run()

if __name__ == '__main__':
	main()
