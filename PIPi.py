#!/usr/bin/env python3

import os, sys, subprocess, argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

pd.options.mode.chained_assignment = None


class Pesti():
		def __init__(self,path,outdir,threads,refs,adapters):
			self.path = path
			self.outdir = outdir
			self.threads = threads
			self.refs = refs
			self.adapters = adapters
		
		def run_fastqc(self,path,outdir):
			print("\n")
			if not os.path.exists(os.path.join(outdir,"fastqc")):
				os.mkdir(os.path.join(outdir,"fastqc"))
			for file in os.listdir(path):
				(base, ext) = os.path.splitext(file)
				if ext in (".fastq", ".fq") and "pass_1" in base:
					# Find the corresponding R2 file
					r2_file = base.replace("pass_1", "pass_2") + ext
					if r2_file not in os.listdir(path):
						continue
					r1_path = os.path.join(path, file)
					r2_path = os.path.join(path, r2_file)
					print(f"Running fastqc on {os.path.join(path,file)}")
					subprocess.run(['fastqc','-o',os.path.join(outdir,"fastqc"), r1_path, r2_path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

		def trimming(self, path, outdir, adapters):
			print("\n")
			if not os.path.exists(os.path.join(outdir,"1_trimming")):
				os.mkdir(os.path.join(outdir,"1_trimming"))
			for file in os.listdir(path):
				(base, ext) = os.path.splitext(file)
				if ext in (".fastq", ".fq") and "pass_1" in base:
					# Find the corresponding R2 file
					r2_file = base.replace("pass_1", "pass_2") + ext
					if r2_file not in os.listdir(path):
						continue
					r1_path = os.path.join(path, file)
					r2_path = os.path.join(path, r2_file)
					print(f"Running bbduk on {os.path.join(path,file)}")
					subprocess.run(['bbduk.sh','-Xmx2g',f'in1={r1_path}',f'in2={r2_path}',f'out1={os.path.join(outdir,"1_trimming/",base+".trim.fastq")}',f'out2={os.path.join(outdir,"1_trimming/",r2_file+".trim.fastq")}',f'stats={os.path.join(outdir,"1_trimming",base+".stats.txt")}',f'ref={adapters}','ktrim=r','k=23','mink=11','hdist=1','qtrim=rl','trimq=10','maq=10','minlen=50','threads=48','tpe','tbo'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

		def map2ref(self, path, outdir, threads,refs):
			print("\n")
			if not os.path.exists(os.path.join(outdir,"2_map2ref")):
				os.mkdir(os.path.join(outdir,"2_map2ref"))
			subprocess.run(['bwa','index','-p','PIPi_index_temp',refs],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
			for file in os.listdir(os.path.join(outdir,"1_trimming")):
				(base, ext) = os.path.splitext(file)
				if ext in (".fastq", ".fq") and "pass_1" in base:
					# Find the corresponding R2 file
					r2_file = base.replace("pass_1", "pass_2") + ext
					if r2_file not in os.listdir(path):
						continue
					idx=base.replace("pass_1","")
					r1_path=os.path.join(outdir,"1_trimming",file)
					r2_path=os.path.join(outdir,"1_trimming",r2_file)
					print(f"Mapping reads from {base} and {r2_file} to {self.refs}")
					sam=subprocess.Popen(['bwa','mem','-t',str(threads),'PIPi_index_temp',r1_path, r2_path],stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
					subprocess.run(['samtools','sort', '-o',os.path.join(outdir,'2_map2ref',base+'.sort.bam'),'-'],stdin=sam.stdout,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
					subprocess.run(['samtools','coverage', '-o',os.path.join(outdir,"2_map2ref",base+".sort.cov"),os.path.join(outdir,'2_map2ref',base+'.sort.bam')],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
					subprocess.run(['samtools','index',os.path.join(outdir,'2_map2ref',base+'.sort.bam')])
					# with open(os.path.join(outdir,'2_map2ref',base+'.sort.bedcov'),'wb') as fh:
					# 	subprocess.run(['samtools','bedcov','data/WG_BVDV.20190405.designed.bed',os.path.join(outdir,'2_map2ref',base+'.sort.bam')],stdout=fh,stderr=subprocess.DEVNULL )
					# print(f"Calculating coverage for {file}\n")


		def output_graph_cov(self,dataframe,outdir,base):
			if not os.path.exists(os.path.join(outdir,"graphs")):
				os.mkdir(os.path.join(outdir,"graphs"))
			subset_data=dataframe[["rname","Percentage"]]
			other = 100 - subset_data['Percentage'].sum()
			df = pd.DataFrame([["Other", other]], columns=["rname","Percentage"])
			percentage_appended=subset_data.append(df)
			percentage_appended.set_index('rname', inplace=True)
			sizes=percentage_appended['Percentage']
			# I'm not sure about this function inside here... but it works, so I'm sticking with it.
			def absolute_value(val):
				a= np.round(val/100.*sizes.sum(), 1)
				return a
			plot=percentage_appended.plot.pie(y='Percentage', figsize=(20,20), labeldistance=0.5, autopct=absolute_value)
			plot.get_legend().remove()
			plt.savefig(os.path.join(outdir,'graphs',base+'.png'))

		def parse_cov(self,outdir):
			for file in os.listdir(os.path.join(outdir,"2_map2ref")):
				(base,ext)=os.path.splitext(file)
				if ext == ".cov":
					input_data=pd.read_csv(os.path.join(outdir,"2_map2ref",file), sep="\t")
					total = input_data['numreads'].sum()
					percentage=input_data["numreads"] / total *100
					input_data["Percentage"]=percentage
					sorted_data=input_data.loc[(input_data['Percentage'] >= 5)].sort_values(by="Percentage", ascending=False).round({'Percentage':2})
					sorted_data.to_csv(os.path.join(outdir,'2_map2ref',base+'cov_results.csv'))
					out_format=sorted_data.rename(columns={'#rname':'rname'})
					self.output_graph_cov(out_format,outdir,base)

		def parse_bedcov(self,outdir):
			for file in os.listdir(os.path.join(outdir,"2_map2ref")):
				(base,ext)=os.path.splitext(file)
				if ext == ".bedcov":
					input_data=pd.read_csv(os.path.join(outdir,"2_map2ref",file),sep="\t", header=None)
					total=input_data[6].sum()
					percentage=input_data[6] / total *100
					input_data["Percentage"]=percentage
					sorted_data=input_data.loc[(input_data['Percentage'] >= 1)].sort_values(by="Percentage", ascending=False).round({'Percentage':2})
					sorted_data=sorted_data.rename(columns={3:'rname'})
					print(sorted_data)
					sorted_data.to_csv(os.path.join(outdir,'2_map2ref',base+'bedcov_results.csv'))
					sorted_data['rname']=sorted_data['rname'].str.replace(r'_1.+',"", regex=True)
					subset_data=sorted_data[['rname','Percentage']]
					subset_data=subset_data.groupby('rname', as_index=False)['Percentage'].sum()
					self.output_graph_cov(subset_data,outdir,base+'_bedcov')

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
			self.trimming(self.path,self.outdir, self.adapters)
			self.map2ref(self.path, self.outdir, self.threads,self.refs)
			self.parse_cov(self.outdir)
			# self.parse_bedcov(self.outdir)
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
	required.add_argument("-a", "--adapters", type=str, required=True, help="Location of the adapter sequences fasta file")
	optional.add_argument("-t", "--threads", type=str, required=False, default=8, help="Number of threads to use. Default 8")
	optional.add_argument("-h", "--help", action="help", help="show this help message and exit")


	args=parser.parse_args()
	path=os.path.normpath(args.path)
	outdir=args.outdir
	refs=args.reference
	adapters=args.adapters
	threads=args.threads

	job=Pesti(path,outdir,threads,refs, adapters)
	job.run()

if __name__ == '__main__':
	main()
