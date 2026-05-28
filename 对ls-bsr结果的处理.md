# LS-BSR结果的处理

## 一、数据准备

参考文献：[Visualization of comparative genomic analyses by BLAST score ratio | BMC Bioinformatics | Full Text](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-2)

### **1.LS-BSR工具的基本原理：**

![pic_1](C:\Users\天一\Desktop\guide_for_freshman\1.4_result_dealing\pics\pic_1.png)

### **2.LS-BSR的结果解读：**

![pic_2](C:\Users\天一\Desktop\guide_for_freshman\1.4_result_dealing\pics\pic_2.png)

某个Straphylococcus的CDS与数据库中所有菌的Complete Genome进行比对的结果

的列标签是细菌的Complete Genome，行标签是这一个菌的CDS

每一个数据也就是所谓的BSR值——==特异度与相似度==

### 3.什么是特异度和相似度

**特异度**：该CDS与非本种菌进行比对后的BSR值，这个值越低，则说明这个基因越特异

**相似度**：该CDS与本种菌进行比对后的BSR值，这个值越高，则说明在本种菌种越有代表性

*<u>我们在筛选基因的过程中，希望这个基因的特异度越低，相似度越高</u>*

## 二、处理代码的解读

### **1.数据分离**

首先给ls-bsr结果分开，分成本种菌和非本种菌两个

也可以在使用ls-bsr的时候，直接分开进行比对，不过放一起更方便一些

```python
def divide(aim_path,include_name,exclude_name,spe_name):	#初步进行筛选，分成基因-本菌和非本菌两个文件
	home_path = os.getcwd()
	df = pd.read_csv(aim_path, sep='\t')
	col_names_list = [col for col in df.columns.tolist() if col != 'Unnamed: 0']    ###【1】
	include_names = []
	exclude_names = []
	for b in col_names_list: 
		if b.split('_')[0] == spe_name:
			include_names.append(b)
		else:
			exclude_names.append(b)
	gene_id = df.iloc[:, 0].tolist()
	inclu_df = df[include_names]#同种类菌比对的结果-相似度
	exclu_df = df[exclude_names]#其余种类菌比对的结果-特异度
	inclu_df.insert(0,'gene_id',gene_id)                                            ###【2】
	exclu_df.insert(0,'gene_id',gene_id)
	include_file_path = os.path.join(home_path,include_name)
	exclude_file_path = os.path.join(home_path,exclude_name)
	inclu_df.to_csv(include_file_path, sep='\t',index = False)
	exclu_df.to_csv(exclude_file_path, sep='\t',index = False)
	print('{} has finished'.format(spe_name))
	return print('\n###Admiral , we have the plans!')
```

【1】因为ls-bsr生成的文件，使用pandas读取列标签的时候会出现一个空缺（Unnamed: 0），所以需要排除此

【2】在生成的新文件exclude.txt和include.txt中，使用这个代码更正了【1】中的小瑕疵

### 2.数据初步整理

这一步我们把文件中的相似度与特异度进行计算，进而获取以下特征以供后续分析

==sim_average==：相似度平均值  ==spe_average==：特异度平均值

==sim_min==：相似度最小值  ==spe_max==：特异度最大值  

==gap==：最小值和最大值的差值(“间距”)  ==comprehensive_score==：根据权重计算得到的分数

```python
home_path = os.getcwd()
exclu_path = os.path.join(home_path,exclude_name)
inclu_path = os.path.join(home_path,include_name)
sim_df = pd.read_csv(inclu_path , sep='\t')
spe_df = pd.read_csv(exclu_path, sep='\t')
#对include数据进行分析
sim_numeric_cols = sim_df.select_dtypes(include=['number']).columns
gene_id = sim_df['gene_id'].tolist()
sim_average = sim_df[sim_numeric_cols].mean(axis=1, skipna=True).fillna(0).tolist() ##计算平均值
sim_min = sim_df[sim_numeric_cols].min(axis=1, skipna=True).fillna(0).tolist()  ##获取最小值
#对exclude数据进行分析
spe_numeric_cols = spe_df.select_dtypes(include=['number']).columns
spe_average = spe_df[spe_numeric_cols].mean(axis=1, skipna=True).fillna(0).tolist() ##计算平均值
spe_max = spe_df[spe_numeric_cols].max(axis=1, skipna=True).fillna(0).tolist() ##获取最大值
gap = [a - b for a, b in zip(sim_min, spe_max)]  #计算
#综合分析后将数据初步整合
pre_df = pd.DataFrame({
	'gene_id':gene_id,
	'sim_average':sim_average,
	'spe_average':spe_average,
	'sim_min':sim_min,
	'spe_max':spe_max,
	'gap':gap
    })
```

![pic_7](C:\Users\天一\Desktop\guide_for_freshman\1.4_result_dealing\pics\pic_7.png)

### ==**3.数据的计算分析**==

#### 3.1 “交叉现象”

![pic3](C:\Users\天一\Desktop\guide_for_freshman\1.4_result_dealing\pics\pic3.png)

​    **交叉：指的是这个基因的特异度最大值大于相似度最小值，进而导致根据其exclude和include文件绘制的折线图出现了重合的部分**

现在市面上存在的核酸试剂盒存在一定的假阳性或者假阴性的情况，这个的原理与这个交叉现象基本吻合。

所以在进行筛选的时候，我们需要尽可能选择交叉没那么严重甚至完全没有交叉的基因

**在实际操作中，我们对基因进行双重确认，需要对应的基因同时满足两个条件，才可以被剔除**

1.该基因的相似度小于特异度最大值的数量超过5%

2.该基因的特异度大于相似度最小值的数量超过5%

```python
spe_total_count = len(spe_numeric_cols.tolist())
sim_total_count = len(sim_numeric_cols.tolist())
sim_min_series = pd.Series(sim_min, index=spe_df.index)
spe_max_series = pd.Series(spe_max, index=sim_df.index)
sim_min_array = sim_min_series.to_numpy()
spe_max_array = spe_max_series.to_numpy()
sim_min_array_2d = sim_min_array.reshape(-1, 1)
spe_max_array_2d = spe_max_array.reshape(-1, 1)
#相似度的处理
sim_row_mask = spe_df[spe_numeric_cols] > sim_min_array_2d
sim_count_less = sim_row_mask.sum(axis=1)
sim_to_drop = sim_count_less >= spe_total_count*0.05
sim_dropped_indices = pre_df[sim_to_drop].index.tolist()
sim_dropped_first_elements = pre_df[sim_to_drop]['gene_id'].tolist()
#特异度的处理
spe_row_mask = sim_df[sim_numeric_cols] < spe_max_array_2d	
spe_count_less = spe_row_mask.sum(axis=1)
spe_to_drop = spe_count_less >= sim_total_count*0.05
spe_dropped_indices = pre_df[spe_to_drop].index.tolist()
spe_dropped_first_elements = pre_df[spe_to_drop]['gene_id'].tolist()
#满足两个条件才可以被剔除
set1 = set(sim_dropped_first_elements)
set2 = set(spe_dropped_first_elements)
intersection_set = set1 & set2
dropped_first_elements = list(intersection_set)
raw_df = pre_df[~pre_df['gene_id'].isin(dropped_first_elements)].copy()
def check_list_empty(input_list):
	return "\n###you are clear!" if not input_list else "\n###Gentlemen one of you betrayed the empire!"
print(check_list_empty(dropped_first_elements))
print('\n')
print(dropped_first_elements)  #输出被剔除的基因
print('\n')
print(raw_df)
```

#### 3.2 min-max标准化

Min-Max 标准化（又称**最小 - 最大归一化 / 离差标准化**）是一种**线性数据缩放方法**，核心是将原始数据按比例映射到 **[0, 1]**（或其他指定区间），==消除量纲与尺度差异，同时保留数据的相对分布关系==。

![pic4](C:\Users\天一\Desktop\guide_for_freshman\1.4_result_dealing\pics\pic4.png)

![pic5](C:\Users\天一\Desktop\guide_for_freshman\1.4_result_dealing\pics\pic5.png)

****

**标准化代码**

```python
#min-max标准化程序
def min_max_scaler(series):
    """Min-Max归一化：(值 - 最小值)/(最大值 - 最小值)"""
    return (series - series.min()) / (series.max() - series.min())
```

```python
#进行标准化操作
raw_df["similarity_scaled"] = min_max_scaler(raw_df["sim_average"])
raw_df["specificity_scaled"] = 1 - min_max_scaler(raw_df["spe_average"])
raw_df['gap_scaled'] = min_max_scaler(raw_df['gap'])
#权重设置——根据具体情况进行选择
weight_similarity = 0.4  
weight_specificity = 0.4
weight_gap = 0.2
raw_df["comprehensive_score"] = (raw_df["similarity_scaled"] * weight_similarity) + \														(raw_df["specificity_scaled"] * weight_specificity) + \
								(raw_df["gap_scaled"] * weight_gap)
df_sorted = raw_df.sort_values(by="comprehensive_score", ascending=False).reset_index(drop=True)
df_sorted.to_csv(
	"gene_data.txt",          # 这个文件名就固定了
    sep='\t',                 # 列分隔符（制表符）
    index=False,              # 不保存行索引（避免多余列）
    encoding='utf-8'          # 避免中文乱码
	)
return print('\n###you may fire when ready.')
```

![pic6](C:\Users\天一\Desktop\guide_for_freshman\1.4_result_dealing\pics\pic6.png)

### 4.散点图的绘制

**具体意义：在参加比对的基因中，展示出选出的具有代表性的三个或一个基因**

![Salmonella](C:\Users\天一\Desktop\guide_for_freshman\1.4_result_dealing\pics\Salmonella.png)

```python
def scatter(spe_name):
	# 设置字体为新罗马字体
	plt.rcParams['font.family'] = 'Times New Roman'
	plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题
	# 读取数据
	data = pd.DataFrame([line.strip().split('\t') for line in open('gene_data.txt', 'r').readlines()])
	data.columns = data.iloc[0]
	data = data[1:].reset_index(drop=True)
	# 转换数据类型
	numeric_cols = ['sim_average', 'spe_average', 'sim_min', 'spe_max', 'similarity_scaled', 'specificity_scaled', 'comprehensive_score']
	for col in numeric_cols:
		data[col] = data[col].astype(float)
	# 提取前三的基因
	top3 = data.head(3)
	first = data.iloc[0:1]
	# 创建画布
	fig, ax = plt.subplots(figsize=(12, 8))
	# 1. 首先绘制所有散点（作为底层，确保不被遮挡）
	scatter = ax.scatter(data['spe_average'], data['sim_average'],
    	                c=data['comprehensive_score'], cmap='viridis',
        	            alpha=0.6, s=50, edgecolors='gray', linewidth=0.5, zorder=1)
	# 2. 绘制前三的基因标记（使用空心样式，确保不遮挡底层点）
	ax.scatter(top3['spe_average'], top3['sim_average'],
           	s=150, facecolors='none', edgecolors='blue', marker='o', linewidth=2,
           	label='Top 3 Genes', zorder=2)
	# 3. 重点标记第一名基因（使用实心标记，放在最上层，不添加文字标注）
	ax.scatter(first['spe_average'], first['sim_average'],
           	s=200, c='none', marker='*', edgecolors='red', linewidth=2,
           	label='Top 1 Gene', zorder=3)
	# 设置坐标轴标签和标题
	ax.set_xlabel('Specificity Average', fontsize=14, fontweight='bold')
	ax.set_ylabel('Similarity Average', fontsize=14, fontweight='bold')
	ax.set_title('Scatter Plot of Gene Similarity and Specificity',
             	fontsize=16, fontweight='bold', pad=20)
	# 添加颜色条
	cbar = plt.colorbar(scatter)
	cbar.set_label('Comprehensive Score', fontsize=12, fontweight='bold')
	# 添加网格和图例
	ax.grid(True, alpha=0.3, linestyle='--')
	ax.legend(loc='upper right', fontsize=12, framealpha=0.9)
	# 调整布局
	plt.tight_layout()
	# 保存图片
	save_name = '{}.png'.format(spe_name)
	plt.savefig(save_name, dpi=300, bbox_inches='tight')
	plt.close()
	#输出有关信息
	print(f"\nTop 1 Gene: {first.iloc[0]['gene_id']}")
	print(f"Top 3 Genes: {', '.join(top3['gene_id'].tolist())}")
	return print("\n###Oh it's beautiful!")
```

### 5.基因比对结果的可视化

![Salmonella_gene1986](C:\Users\天一\Desktop\guide_for_freshman\1.4_result_dealing\pics\Salmonella_gene1986.png)

这个图的左边显示的是非本种菌，右边显示的是本种菌，纵坐标表示的是细菌数量，横坐标表示的是特异度/相似度

```python
def draw(include_name,exclude_name,draw_2,spe_name):
	home_path = os.getcwd()
	os.mkdir(draw_2)
	draw_2 = os.path.join(home_path,draw_2)
	exclu_path = os.path.join(home_path,exclude_name)
	inclu_path = os.path.join(home_path,include_name)
	df_1 = pd.read_csv(inclu_path, sep='\t')
	df_2 = pd.read_csv(exclu_path, sep='\t')
	data_1 = df_1.apply(pd.to_numeric, errors='coerce')
	data_2 = df_2.apply(pd.to_numeric, errors='coerce')
	gene_df = pd.read_csv("gene_data.txt", sep="\t")
	# 按综合得分降序排序，取前三基因
	top3_genes = gene_df.nlargest(3, "comprehensive_score")["gene_id"].tolist()
	for idx,val in enumerate(top3_genes):
		save_name = '{}.png'.format(val)   #基因名字命名文件
		save_path = os.path.join(draw_2,save_name)
		in_lst = df_1[df_1['gene_id'] == val].values.tolist()[0][1:]
		ex_lst = df_2[df_2['gene_id'] == val].values.tolist()[0][1:]
		ex_lst.insert(0, -0.005)
		counts_in = Counter(in_lst)
		counts_ex = Counter(ex_lst)
		# 提取值和对应的计数
		values_in = sorted(counts_in.keys())
		counts_list_in = [counts_in[value] for value in values_in]
		values_ex = sorted(counts_ex.keys())
		counts_list_ex = [counts_ex[value] for value in values_ex]
		counts_list_ex[0] = 0
		fig,ax1 = plt.subplots()
		plt.rcParams['font.family'] = 'Times New Roman'
		plt.rcParams['figure.dpi'] = 500
		# 绘制折线图
		ax1.plot(values_ex,counts_list_ex,color='#F6631C', linestyle='-',linewidth=2)
		ax1.set_ylabel('amount of other bacteria',color='#F6631C',fontsize=25)
		ax1.tick_params(axis='y',colors='#F6631C',width = 2)
		ax1.minorticks_on()
		ax1.yaxis.set_minor_locator(MultipleLocator(2500))  # Y轴次要刻度间隔
		ax1.xaxis.set_minor_locator(MultipleLocator(0.05))  # X轴次要刻度间隔
		ax1.tick_params(which='minor', axis='both', length=0)  # 设置次要刻度线长度为0
		ax2 = ax1.twinx()
		ax2.plot(values_in,counts_list_in,color='#6D65A3', linestyle='-', linewidth=2)
		ax2.set_ylabel('amount of {}'.format(spe_name),color='#6D65A3',fontsize=25)
		ax2.tick_params(axis='y',colors='#6D65A3',width = 2)
		plt.title(val,fontsize=25,fontweight='bold')
		ax = plt.gca()
		plt.xlim(-0.025, 1.025)
		ax.set_xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
		plt.xticks(fontsize=12, fontweight="bold")
		for tick in ax1.get_xticklabels():
			tick.set_weight("bold")
			tick.set_size(12)
		for tick in ax1.get_yticklabels():
			tick.set_weight("bold")
			tick.set_size(12)
		for tick in ax2.get_yticklabels():
			tick.set_weight("bold")
			tick.set_size(12)
		# 显示网格线
		ax1.set_facecolor('#F8F8F8')
		ax1.grid(which='major', axis='both',color='white',linestyle='-',linewidth=2.5,alpha=1)
		ax1.grid(which='minor', axis='both', linestyle='-',color='white', alpha=1, linewidth=1.5)
		# 显示图形
		plt.tight_layout()
		plt.savefig(save_path)
		plt.cla()
		plt.close("all")
	return print('\n###Rebellions are built on hope!')
```































































