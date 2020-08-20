# -*- coding: utf-8 -*-

import argparse
import traceback
import logging
import math
import numbers
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

from .biostats import ppoints

def manhattan(assoc,out=False,gap=10,ax=False,cmap=False,cmap_var=2,
			  col_chr="CHR",col_bp="BP",col_p='P',col_snp="SNP",show=False,
			  title=False,xtick_size=10,ytick_size=10,xrotation=0,yrotation=0,label_size=15,title_size=20,
			  suggestiveline=-np.log10(1e-5), genomewideline=-np.log10(5e-8),**kwargs):

	if not (isinstance(cmap_var, numbers.Number) or isinstance(cmap_var, list)):
		raise Exception("[ERROR]: cmap_var should either be list or number.")

	list_color = list()
	if cmap:
		try:
			if isinstance(cmap_var, numbers.Number):
				step = int(len(cmap.colors) / cmap_var)

				for i, color in enumerate(cmap.colors):
					if i % step == 0:
						list_color.append(colors.to_hex(color))
			else:
				list_color = cmap_var

		except AttributeError:
			if isinstance(cmap_var, numbers.Number):
				step = int(256 / cmap_var)

				for i in range(256):
					if i % step == 0:
						list_color.append(colors.to_hex(cmap(i)[:3]))
			else:
				list_color = cmap_var
	else:
		cmap = plt.get_cmap("Greys_r")
		if isinstance(cmap_var, numbers.Number):
			if cmap_var == 2:
				list_color = [colors.to_hex(cmap(0)[:3]), colors.to_hex(cmap(80)[:3])]
			else:
				step = int(256 / cmap_var)

				for i in range(256):
					if i % step == 0:
						list_color.append(colors.to_hex(cmap(i)[:3]))
		else:
			list_color = cmap_var

	if not (ax or show or out):
		raise Exception("[ERROR]: Either of the ax, show, and out must have a value.")
	isAx = not ax
	if isinstance(assoc, str):
		df_assoc = pd.read_csv(assoc, header=0, delim_whitespace=True)
	elif isinstance(assoc, pd.DataFrame):
		df_assoc = assoc
	else:
		raise Exception("[ERROR]: assoc must be either string(path) or pandas.DataFrame.")
	
	if col_chr not in df_assoc.columns:	raise Exception("[ERROR]: Column '{0}' not found!".format(col_chr))
	if col_bp not in df_assoc.columns:	raise Exception("[ERROR]: Column '{0}' not found!".format(col_bp))
	if col_p not in df_assoc.columns:	raise Exception("[ERROR]: Column '{0}' not found!".format(col_p))
	if col_snp not in df_assoc.columns:	print("[WARNING]: Column '{0}' not found!".format(col_snp))

	df_assoc[col_chr] = df_assoc[col_chr].astype('category')
	df_assoc[col_bp] = df_assoc[col_bp].astype(int)
	df_assoc[col_p] = df_assoc[col_p].astype(float)
	df_assoc = df_assoc.sort_values([col_chr,col_bp])

	chr_len = list()
	for cChr in df_assoc[col_chr].unique():
		chr_len.append(len(df_assoc[df_assoc[col_chr]==cChr]))
	weight_gap = int(min(chr_len)/100)
	
	list_ind = list()
	for cChr in df_assoc[col_chr].unique():
		if len(list_ind) == 0:
			last_ind = 0
		else:
			last_ind = (list_ind[-1]+1)+(gap*weight_gap)
		
		list_ind += [last_ind+num for num in range(len(df_assoc[df_assoc[col_chr]==cChr]))]
		
	df_assoc['IND'] = list_ind
	df_assoc["LOG_P"] = -np.log10(df_assoc[col_p])
	
	if isAx:
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (12,12))
	else:
		fig = ax.figure
	
	x_ticks,x_labels = list(),list()


	for i, cChr in enumerate(df_assoc[col_chr].unique()):
		ind = df_assoc[df_assoc[col_chr]==cChr]['IND']
		log_p = df_assoc[df_assoc[col_chr]==cChr]['LOG_P']

		ax.scatter(ind, log_p, marker='.',s=5,color=list_color[i%cmap_var],**kwargs)
		x_ticks.append(ind.iloc[0]+(ind.iloc[-1]-ind.iloc[0])/2)
		x_labels.append(cChr)
	
	x_padding = len(list_ind)/20
	xlim_min = list_ind[0]-x_padding
	xlim_max = list_ind[-1]+x_padding
	
	if suggestiveline:
		ax.plot([xlim_min,xlim_max],[suggestiveline,suggestiveline],'b-')
	if genomewideline:
		ax.plot([xlim_min,xlim_max],[genomewideline,genomewideline],'r-')
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)

	ax.set_xticks(x_ticks)
	ax.set_xticklabels(x_labels)
	# ax.set_ylim(bottom=0,top=math.floor(1+max(df_assoc["LOG_P"])+(max(df_assoc["LOG_P"])/100)))
	# ax.set_ylim(bottom=0,top=max(df_assoc["LOG_P"]))
	ax.set_ylim(bottom=0)
	ax.set_xlim([xlim_min,xlim_max])

	ax.tick_params(axis='x', labelsize=xtick_size, labelrotation=xrotation)
	ax.tick_params(axis='y', labelsize=ytick_size, labelrotation=yrotation)
	
	ax.set_xlabel("Chromosomes",fontsize=label_size)
	ax.set_ylabel(r'$-log_{10}(p)$',fontsize=label_size)

	if title:
		ax.set_title(title,fontsize=title_size)
	
	if isAx:
		fig.tight_layout()

	if show:
		plt.show()
		
	if out:
		plt.savefig(out,format="png")
		
	if isAx:
		plt.clf()
		plt.close()

def qqplot(assoc,out=False,col_p='P',show=False,ax=False,
		   title=False,xtick_size=10,ytick_size=10,xrotation=0,yrotation=0,label_size=15,title_size=20,
		   **kwargs):
	if not (ax or show or out):
		raise Exception("[ERROR]: Either of the ax, show, and out must have a value.")
		
	isAx = not ax
	p_vals = None
	if isinstance(assoc, str):
		df_assoc = pd.read_csv(assoc, header=0, delim_whitespace=True, dtype={1:'category', 2:int, 4:float, 5:float, 7:float, 8:float, 9:float})
		p_vals = df_assoc[col_p].dropna()
		p_vals = p_vals[(0<p_vals)&(p_vals<1)]
		
	elif isinstance(assoc, pd.DataFrame):
		p_vals = assoc[col_p].dropna()
		p_vals = p_vals[(0<p_vals)&(p_vals<1)]
	elif isinstance(assoc,pd.Series):
		p_vals = assoc.dropna()
		p_vals = p_vals[(0<p_vals)&(p_vals<1)]
	else:
		p_vals = [ele for ele in assoc if (0<ele) and (ele<1)]
	
	observed = -np.log10(np.sort(np.array(p_vals)))
	expected = -np.log10(ppoints(len(p_vals)))
	
	x_padding = (np.nanmax(expected)-np.nanmin(expected))/12
	y_padding = (np.nanmax(observed)-np.nanmin(observed))/12
	
	
	if isAx:
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (12,12))
	else:
		fig = ax.figure

	ax.scatter(expected,observed,c='k',**kwargs)
	
	xlim_min = np.nanmin(expected)-x_padding
	xlim_max = np.nanmax(expected)+x_padding
	ylim_min = np.nanmin(observed)-y_padding
	ylim_max = np.nanmax(observed)+y_padding
	
	max_lim = xlim_max if xlim_max<ylim_max else ylim_max
	min_lim = xlim_min if xlim_min>ylim_min else ylim_min
	ax.plot([min_lim,max_lim],[min_lim,max_lim],'r-')
	
	ax.set_xlim([xlim_min, xlim_max])
	ax.set_ylim([ylim_min, ylim_max])
	ax.set_xlabel("Expected $-log_{10}(p)$",fontsize=label_size)
	ax.set_ylabel("Observed $-log_{10}(p)$",fontsize=label_size)

	ax.tick_params(axis='x', labelsize=xtick_size, labelrotation=xrotation)
	ax.tick_params(axis='y', labelsize=ytick_size, labelrotation=yrotation)

	if title:
		ax.set_title(title,fontsize=title_size)

	if isAx:
		fig.tight_layout()

	if show:
		plt.show()
		
	if out:
		plt.savefig(out,format="png")
		
	if isAx:
		plt.clf()
		plt.close()
	
if __name__ == "__main__":
	FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
	logging.basicConfig(level=logging.INFO, format=FORMAT)
	logger = logging.getLogger(__name__)
	
	parser = argparse.ArgumentParser(description='PGMplatform')
	
	parser.add_argument('--assoc', type=str,help=' ',default='./',required=True)
	parser.add_argument('--out', type=str,help=' ',default=None,required=False)
	parser.add_argument('--gap', type=int,help=' ',default=10,required=False)
	parser.add_argument('--plot', type=str,help='[manhattan/qqplot]',required=True)
	parser.add_argument('--show', action='store_true',required=False)
	
	args = parser.parse_args()
	
	logger.info(args)

	try:
		if args.plot in ["Manhattan", "manhattan"]:
			if args.out == None:
				manhattan(args.assoc,"./Manhattan.png",args.gap,show=args.show)
			else:
				manhattan(args.assoc,args.out,args.gap,show=args.show)
		elif args.plot in ["QQplot","qqplot","qq","QQ","QQPlot"]:
			if args.out == None:
				qqplot(args.assoc,"./QQplot.png")
			else:
				qqplot(args.assoc,args.out)
					 
	except Exception:
		logger.error(traceback.format_exc())
		logger.error("qqman.py failure on arguments: {}".format(args))