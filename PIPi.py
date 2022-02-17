#!/usr/bin/env python3

import os, sys, subprocess, argparse, csv
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from distutils.spawn import find_executable

pd.options.mode.chained_assignment = None


class Pesti():
		def __init__(self,path,outdir,threads,refs,adapters,bed):
			self.path = path
			self.outdir = outdir
			self.threads = threads
			self.refs = refs
			self.adapters = adapters
			self.bed=bed


		def run_fastqc(self,path,outdir):
			print("\n")
			if not os.path.exists(os.path.join(outdir,"fastqc")):
				os.mkdir(os.path.join(outdir,"fastqc"))
			for file in os.listdir(path):
				(base,ext)=os.path.splitext(file)
				if ext == ".fastq":
					print(f"Running fastqc on {os.path.join(path,file)}")
					subprocess.run(['fastqc','-o',os.path.join(outdir,"fastqc"), os.path.join(path,file)],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

		def trimming(self, path, outdir, adapters):
			print("\n")
			if not os.path.exists(os.path.join(outdir,"1_trimming")):
				os.mkdir(os.path.join(outdir,"1_trimming"))
			for file in os.listdir(path):
				(base,ext)=os.path.splitext(file)
				if ext == ".fastq":
					print(f"Running bbduk on {os.path.join(path,file)}")
					subprocess.run(['bbduk.sh','-Xmx2g',f'in1={os.path.join(path,file)}',f'out1={os.path.join(outdir,"1_trimming/",base+".trim.fastq")}',f'stats={os.path.join(outdir,"1_trimming",base+".stats.txt")}',f'ref={adapters}','ktrim=r','k=23','mink=11','hdist=1','qtrim=rl','trimq=10','maq=10','minlen=50','threads=48','tpe','tbo'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

		def map2ref(self, path, outdir, threads,refs,bed):
			print("\n")
			if not os.path.exists(os.path.join(outdir,"2_map2ref")):
				os.mkdir(os.path.join(outdir,"2_map2ref"))
			subprocess.run(['bwa','index','-p','PIPi_index_temp',refs],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
			for file in os.listdir(os.path.join(outdir,"1_trimming")):
				(base,ext)=os.path.splitext(file)
				if ext == ".fastq":
					print(f"Mapping reads from {file} to {refs}")
					sam=subprocess.Popen(['bwa','mem','-t',str(threads),'PIPi_index_temp',os.path.join(outdir,"1_trimming",file)],stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
					subprocess.run(['samtools','sort', '-o',os.path.join(outdir,'2_map2ref',base+'.sort.bam'),'-'],stdin=sam.stdout,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
					print(f"Calculating coverage for {file}\n")
					subprocess.run(['samtools','coverage', '-o',os.path.join(outdir,"2_map2ref",base+".sort.cov"),os.path.join(outdir,'2_map2ref',base+'.sort.bam')],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
					subprocess.run(['samtools','index',os.path.join(outdir,'2_map2ref',base+'.sort.bam')])
					with open(os.path.join(outdir,'2_map2ref',base+'.sort.bedcov'),'w') as fh:
						subprocess.run(['samtools','bedcov',bed,os.path.join(outdir,'2_map2ref',base+'.sort.bam')],stdout=fh)

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
					total=input_data[4].sum()
					percentage=input_data[4] / total *100
					input_data["Percentage"]=percentage
					sorted_data=input_data.loc[(input_data['Percentage'] >= 1)].sort_values(by="Percentage", ascending=False).round({'Percentage':2})
					sorted_data=sorted_data.rename(columns={3:'rname'})
					# print(sorted_data)
					sorted_data.to_csv(os.path.join(outdir,'2_map2ref',base+'bedcov_results.csv'))
					sorted_data['rname']=sorted_data['rname'].str.replace(r'_1.+',"", regex=True)
					subset_data=sorted_data[['rname','Percentage']]
					subset_data=subset_data.groupby('rname', as_index=False)['Percentage'].sum()
					self.output_graph_cov(subset_data,outdir,base+'_bedcov')

		def table_summarise(self,outdir):
			if not os.path.exists(os.path.join(outdir,"summary")):
				os.mkdir(os.path.join(outdir,"summary"))
			top=[['Index','Name', 'Startpos', 'endpos','numreads','covbases','coverage','meandepth','meanbaseq','meannmapq','Percentage','Sample']]
			all=[['Index','Name', 'Startpos', 'endpos','numreads','covbases','coverage','meandepth','meanbaseq','meannmapq','Percentage','Sample']]
			for file in os.listdir(os.path.join(outdir,"2_map2ref")):
				if file[-19:] == "sortcov_results.csv":
					covdf=pd.read_csv(os.path.join(outdir,"2_map2ref",file))
					sort_covdf=covdf.sort_values(by="Percentage", ascending=False)
					for k,row in sort_covdf.iterrows():
						if k == 0:
							ls=sort_covdf.loc[k].values.flatten().tolist()
							ls.append(file[:-25])
							top.append(ls)
							all.append(ls)
						else:
							ls=sort_covdf.loc[k].values.flatten().tolist()
							ls.append(file[:-25])
							all.append(ls)
			with open(os.path.join(outdir,"summary","Top_hits_per_sample.csv"), "w") as f:
				writer = csv.writer(f)
				writer.writerows(top)
			with open(os.path.join(outdir,"summary","All_hits_per_sample.csv"), "w") as f:
				writer = csv.writer(f)
				writer.writerows(all)

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
			self.map2ref(self.path, self.outdir, self.threads,self.refs,self.bed)
			self.parse_cov(self.outdir)
			self.parse_bedcov(self.outdir)
			self.table_summarise(self.outdir)
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
	required.add_argument("-b", "--bed_file", type=str, required=True, help="Location of the bed file of target regions")
	optional.add_argument("-t", "--threads", type=str, required=False, default=8, help="Number of threads to use. Default 8")
	optional.add_argument("-h", "--help", action="help", help="show this help message and exit")


	args=parser.parse_args()
	path=os.path.normpath(args.path)
	outdir=args.outdir
	refs=args.reference
	adapters=args.adapters
	bed=args.bed_file
	threads=args.threads

	job=Pesti(path,outdir,threads,refs,adapters,bed)
	job.run()

if __name__ == '__main__':
	main()
