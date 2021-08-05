cola Report for Hierarchical Partitioning - 'Golub_leukemia'
==================

**Date**: 2021-07-22 16:02:47 CEST, **cola version**: 1.9.4

----------------------------------------------------------------

<style type='text/css'>

body, td, th {
   font-family: Arial,Helvetica,sans-serif;
   background-color: white;
   font-size: 13px;
  max-width: 800px;
  margin: auto;
  margin-left:210px;
  padding: 0px 10px 0px 10px;
  border-left: 1px solid #EEEEEE;
  line-height: 150%;
}

tt, code, pre {
   font-family: 'DejaVu Sans Mono', 'Droid Sans Mono', 'Lucida Console', Consolas, Monaco, 

monospace;
}

h1 {
   font-size:2.2em;
}

h2 {
   font-size:1.8em;
}

h3 {
   font-size:1.4em;
}

h4 {
   font-size:1.0em;
}

h5 {
   font-size:0.9em;
}

h6 {
   font-size:0.8em;
}

a {
  text-decoration: none;
  color: #0366d6;
}

a:hover {
  text-decoration: underline;
}

a:visited {
   color: #0366d6;
}

pre, img {
  max-width: 100%;
}
pre {
  overflow-x: auto;
}
pre code {
   display: block; padding: 0.5em;
}

code {
  font-size: 92%;
  border: 1px solid #ccc;
}

code[class] {
  background-color: #F8F8F8;
}

table, td, th {
  border: 1px solid #ccc;
}

blockquote {
   color:#666666;
   margin:0;
   padding-left: 1em;
   border-left: 0.5em #EEE solid;
}

hr {
   height: 0px;
   border-bottom: none;
   border-top-width: thin;
   border-top-style: dotted;
   border-top-color: #999999;
}

@media print {
   * {
      background: transparent !important;
      color: black !important;
      filter:none !important;
      -ms-filter: none !important;
   }

   body {
      font-size:12pt;
      max-width:100%;
   }

   a, a:visited {
      text-decoration: underline;
   }

   hr {
      visibility: hidden;
      page-break-before: always;
   }

   pre, blockquote {
      padding-right: 1em;
      page-break-inside: avoid;
   }

   tr, img {
      page-break-inside: avoid;
   }

   img {
      max-width: 100% !important;
   }

   @page :left {
      margin: 15mm 20mm 15mm 10mm;
   }

   @page :right {
      margin: 15mm 10mm 15mm 20mm;
   }

   p, h2, h3 {
      orphans: 3; widows: 3;
   }

   h2, h3 {
      page-break-after: avoid;
   }
}
</style>




## Summary



First the variable is renamed to `res_rh`.


```r
res_rh = rh
```



The partition hierarchy and all available functions which can be applied to `res_rh` object.


```r
res_rh
```

```
#> A 'HierarchicalPartition' object with 'ATC:skmeans' method.
#>   On a matrix with 4116 rows and 72 columns.
#>   Performed in total 1350 partitions.
#>   There are 6 groups under the following parameters:
#>     - min_samples: 6
#>     - mean_silhouette_cutoff: 0.9
#>     - min_n_signatures: 102 (signatures are selected based on:)
#>       - fdr_cutoff: 0.05
#>       - group_diff (scaled values): 0.5
#> 
#> Hierarchy of the partition:
#>   0, 72 cols
#>   |-- 01, 35 cols, 650 signatures
#>   |   |-- 011, 11 cols (b)
#>   |   |-- 012, 13 cols, 4 signatures (c)
#>   |   `-- 013, 11 cols (b)
#>   |-- 02, 24 cols, 134 signatures
#>   |   |-- 021, 13 cols, 35 signatures (c)
#>   |   `-- 022, 11 cols (b)
#>   `-- 03, 13 cols, 6 signatures (c)
#> Stop reason:
#>   b) Subgroup had too few columns.
#>   c) There were too few signatures.
#> 
#> Following methods can be applied to this 'HierarchicalPartition' object:
#>  [1] "all_leaves"            "all_nodes"             "cola_report"           "collect_classes"      
#>  [5] "colnames"              "compare_signatures"    "dimension_reduction"   "functional_enrichment"
#>  [9] "get_anno_col"          "get_anno"              "get_children_nodes"    "get_classes"          
#> [13] "get_matrix"            "get_signatures"        "is_leaf_node"          "max_depth"            
#> [17] "merge_node"            "ncol"                  "node_info"             "node_level"           
#> [21] "nrow"                  "rownames"              "show"                  "split_node"           
#> [25] "suggest_best_k"        "test_to_known_factors" "top_rows_heatmap"      "top_rows_overlap"     
#> 
#> You can get result for a single node by e.g. object["01"]
```

The call of `hierarchical_partition()` was:


```
#> hierarchical_partition(data = m, anno = anno[, c("ALL.AML"), drop = FALSE], anno_col = c(ALL = "red", 
#>     AML = "blue"), cores = 4)
```

Dimension of the input matrix:


```r
mat = get_matrix(res_rh)
dim(mat)
```

```
#> [1] 4116   72
```

All the methods that were tried:


```r
res_rh@param$combination_method
```

```
#> [[1]]
#> [1] "ATC"     "skmeans"
```

### Density distribution

The density distribution for each sample is visualized as one column in the following heatmap.
The clustering is based on the distance which is the Kolmogorov-Smirnov statistic between two distributions.




```r
library(ComplexHeatmap)
densityHeatmap(mat, top_annotation = HeatmapAnnotation(df = get_anno(res_rh), 
    col = get_anno_col(res_rh)), ylab = "value", cluster_columns = TRUE, show_column_names = FALSE,
    mc.cores = 1)
```

![plot of chunk density-heatmap](figure_cola/density-heatmap-1.png)



Some values about the hierarchy:


```r
all_nodes(res_rh)
```

```
#> [1] "0"   "01"  "011" "012" "013" "02"  "021" "022" "03"
```

```r
all_leaves(res_rh)
```

```
#> [1] "011" "012" "013" "021" "022" "03"
```

```r
node_info(res_rh)
```

```
#>    id best_method depth best_k n_columns n_signatures p_signatures is_leaf
#> 1   0 ATC:skmeans     1      3        72         2054     0.499028   FALSE
#> 2  01 ATC:skmeans     2      3        35          650     0.157920   FALSE
#> 3 011 not applied     3     NA        11           NA           NA    TRUE
#> 4 012 ATC:skmeans     3      3        13            4     0.000972    TRUE
#> 5 013 not applied     3     NA        11           NA           NA    TRUE
#> 6  02 ATC:skmeans     2      2        24          134     0.032556   FALSE
#> 7 021 ATC:skmeans     3      2        13           35     0.008503    TRUE
#> 8 022 not applied     3     NA        11           NA           NA    TRUE
#> 9  03 ATC:skmeans     2      2        13            6     0.001458    TRUE
```

In the output from `node_info()`, there are the following columns:

- `id`: The node id.
- `best_method`: The best method selected.
- `depth`: Depth of the node in the hierarchy.
- `best_k`: Best number of groups of the partition on that node.
- `n_columns`: Number of columns in the submatrix.
- `n_signatures`: Number of signatures with the `best_k`.
- `p_signatures`: Proportion of hte signatures in total number of rows in the matrix.
- `is_leaf`: Whether the node is a leaf.

Labels of nodes are encoded in a special way. The number of digits
correspond to the depth of the node in the hierarchy and the value of the
digits correspond to the index of the subgroup in the current node, E.g. a label
of “012” means the node is the second subgroup of the partition which is the
first subgroup of the root node.

### Suggest the best k



Following table shows the best `k` (number of partitions) for each node in the
partition hierarchy. Clicking on the node name in the table goes to the
corresponding section for the partitioning on that node.

[The cola vignette](https://jokergoo.github.io/cola_vignettes/cola.html#toc_13)
explains the definition of the metrics used for determining the best
number of partitions.



```r
suggest_best_k(res_rh)
```


|Node              |Best method                                         |Is leaf   |Best k |1-PAC |Mean silhouette |Concordance | #samples|   |
|:-----------------|:---------------------------------------------------|:---------|:------|:-----|:---------------|:-----------|--------:|:--|
|[Node0](#Node0)   |ATC:skmeans                                         |          |3      |1.00  |0.97            |0.98        |       72|** |
|[Node01](#Node01) |ATC:skmeans                                         |          |3      |1.00  |0.95            |0.98        |       35|** |
|Node011-leaf      |<span style='color:grey;'><i>not applied</i></span> |✓ (b)     |       |      |                |            |       11|   |
|Node012-leaf      |ATC:skmeans                                         |✓ (&#99;) |3      |1.00  |0.99            |1.00        |       13|** |
|Node013-leaf      |<span style='color:grey;'><i>not applied</i></span> |✓ (b)     |       |      |                |            |       11|   |
|[Node02](#Node02) |ATC:skmeans                                         |          |2      |1.00  |1.00            |1.00        |       24|** |
|Node021-leaf      |ATC:skmeans                                         |✓ (&#99;) |2      |0.85  |0.91            |0.97        |       13|   |
|Node022-leaf      |<span style='color:grey;'><i>not applied</i></span> |✓ (b)     |       |      |                |            |       11|   |
|Node03-leaf       |ATC:skmeans                                         |✓ (&#99;) |2      |0.71  |0.91            |0.96        |       13|   |


Stop reason: b) Subgroup had too few columns. c) There were too few signatures. 

\*\*: 1-PAC > 0.95, \*: 1-PAC > 0.9


### Partition hierarchy

The nodes of the hierarchy can be merged by setting the `merge_node` parameters. Here we 
control the hierarchy with the `min_n_signatures` parameter. The value of `min_n_signatures` is
from `node_info()`.





<style type='text/css'>



.ui-helper-hidden {
	display: none;
}
.ui-helper-hidden-accessible {
	border: 0;
	clip: rect(0 0 0 0);
	height: 1px;
	margin: -1px;
	overflow: hidden;
	padding: 0;
	position: absolute;
	width: 1px;
}
.ui-helper-reset {
	margin: 0;
	padding: 0;
	border: 0;
	outline: 0;
	line-height: 1.3;
	text-decoration: none;
	font-size: 100%;
	list-style: none;
}
.ui-helper-clearfix:before,
.ui-helper-clearfix:after {
	content: "";
	display: table;
	border-collapse: collapse;
}
.ui-helper-clearfix:after {
	clear: both;
}
.ui-helper-zfix {
	width: 100%;
	height: 100%;
	top: 0;
	left: 0;
	position: absolute;
	opacity: 0;
	filter:Alpha(Opacity=0); 
}

.ui-front {
	z-index: 100;
}



.ui-state-disabled {
	cursor: default !important;
	pointer-events: none;
}



.ui-icon {
	display: inline-block;
	vertical-align: middle;
	margin-top: -.25em;
	position: relative;
	text-indent: -99999px;
	overflow: hidden;
	background-repeat: no-repeat;
}

.ui-widget-icon-block {
	left: 50%;
	margin-left: -8px;
	display: block;
}




.ui-widget-overlay {
	position: fixed;
	top: 0;
	left: 0;
	width: 100%;
	height: 100%;
}
.ui-accordion .ui-accordion-header {
	display: block;
	cursor: pointer;
	position: relative;
	margin: 2px 0 0 0;
	padding: .5em .5em .5em .7em;
	font-size: 100%;
}
.ui-accordion .ui-accordion-content {
	padding: 1em 2.2em;
	border-top: 0;
	overflow: auto;
}
.ui-autocomplete {
	position: absolute;
	top: 0;
	left: 0;
	cursor: default;
}
.ui-menu {
	list-style: none;
	padding: 0;
	margin: 0;
	display: block;
	outline: 0;
}
.ui-menu .ui-menu {
	position: absolute;
}
.ui-menu .ui-menu-item {
	margin: 0;
	cursor: pointer;
	
	list-style-image: url("data:image/gif;base64,R0lGODlhAQABAIAAAAAAAP///yH5BAEAAAAALAAAAAABAAEAAAIBRAA7");
}
.ui-menu .ui-menu-item-wrapper {
	position: relative;
	padding: 3px 1em 3px .4em;
}
.ui-menu .ui-menu-divider {
	margin: 5px 0;
	height: 0;
	font-size: 0;
	line-height: 0;
	border-width: 1px 0 0 0;
}
.ui-menu .ui-state-focus,
.ui-menu .ui-state-active {
	margin: -1px;
}


.ui-menu-icons {
	position: relative;
}
.ui-menu-icons .ui-menu-item-wrapper {
	padding-left: 2em;
}


.ui-menu .ui-icon {
	position: absolute;
	top: 0;
	bottom: 0;
	left: .2em;
	margin: auto 0;
}


.ui-menu .ui-menu-icon {
	left: auto;
	right: 0;
}
.ui-button {
	padding: .4em 1em;
	display: inline-block;
	position: relative;
	line-height: normal;
	margin-right: .1em;
	cursor: pointer;
	vertical-align: middle;
	text-align: center;
	-webkit-user-select: none;
	-moz-user-select: none;
	-ms-user-select: none;
	user-select: none;

	
	overflow: visible;
}

.ui-button,
.ui-button:link,
.ui-button:visited,
.ui-button:hover,
.ui-button:active {
	text-decoration: none;
}


.ui-button-icon-only {
	width: 2em;
	box-sizing: border-box;
	text-indent: -9999px;
	white-space: nowrap;
}


input.ui-button.ui-button-icon-only {
	text-indent: 0;
}


.ui-button-icon-only .ui-icon {
	position: absolute;
	top: 50%;
	left: 50%;
	margin-top: -8px;
	margin-left: -8px;
}

.ui-button.ui-icon-notext .ui-icon {
	padding: 0;
	width: 2.1em;
	height: 2.1em;
	text-indent: -9999px;
	white-space: nowrap;

}

input.ui-button.ui-icon-notext .ui-icon {
	width: auto;
	height: auto;
	text-indent: 0;
	white-space: normal;
	padding: .4em 1em;
}



input.ui-button::-moz-focus-inner,
button.ui-button::-moz-focus-inner {
	border: 0;
	padding: 0;
}
.ui-controlgroup {
	vertical-align: middle;
	display: inline-block;
}
.ui-controlgroup > .ui-controlgroup-item {
	float: left;
	margin-left: 0;
	margin-right: 0;
}
.ui-controlgroup > .ui-controlgroup-item:focus,
.ui-controlgroup > .ui-controlgroup-item.ui-visual-focus {
	z-index: 9999;
}
.ui-controlgroup-vertical > .ui-controlgroup-item {
	display: block;
	float: none;
	width: 100%;
	margin-top: 0;
	margin-bottom: 0;
	text-align: left;
}
.ui-controlgroup-vertical .ui-controlgroup-item {
	box-sizing: border-box;
}
.ui-controlgroup .ui-controlgroup-label {
	padding: .4em 1em;
}
.ui-controlgroup .ui-controlgroup-label span {
	font-size: 80%;
}
.ui-controlgroup-horizontal .ui-controlgroup-label + .ui-controlgroup-item {
	border-left: none;
}
.ui-controlgroup-vertical .ui-controlgroup-label + .ui-controlgroup-item {
	border-top: none;
}
.ui-controlgroup-horizontal .ui-controlgroup-label.ui-widget-content {
	border-right: none;
}
.ui-controlgroup-vertical .ui-controlgroup-label.ui-widget-content {
	border-bottom: none;
}


.ui-controlgroup-vertical .ui-spinner-input {

	
	width: 75%;
	width: calc( 100% - 2.4em );
}
.ui-controlgroup-vertical .ui-spinner .ui-spinner-up {
	border-top-style: solid;
}

.ui-checkboxradio-label .ui-icon-background {
	box-shadow: inset 1px 1px 1px #ccc;
	border-radius: .12em;
	border: none;
}
.ui-checkboxradio-radio-label .ui-icon-background {
	width: 16px;
	height: 16px;
	border-radius: 1em;
	overflow: visible;
	border: none;
}
.ui-checkboxradio-radio-label.ui-checkboxradio-checked .ui-icon,
.ui-checkboxradio-radio-label.ui-checkboxradio-checked:hover .ui-icon {
	background-image: none;
	width: 8px;
	height: 8px;
	border-width: 4px;
	border-style: solid;
}
.ui-checkboxradio-disabled {
	pointer-events: none;
}
.ui-datepicker {
	width: 17em;
	padding: .2em .2em 0;
	display: none;
}
.ui-datepicker .ui-datepicker-header {
	position: relative;
	padding: .2em 0;
}
.ui-datepicker .ui-datepicker-prev,
.ui-datepicker .ui-datepicker-next {
	position: absolute;
	top: 2px;
	width: 1.8em;
	height: 1.8em;
}
.ui-datepicker .ui-datepicker-prev-hover,
.ui-datepicker .ui-datepicker-next-hover {
	top: 1px;
}
.ui-datepicker .ui-datepicker-prev {
	left: 2px;
}
.ui-datepicker .ui-datepicker-next {
	right: 2px;
}
.ui-datepicker .ui-datepicker-prev-hover {
	left: 1px;
}
.ui-datepicker .ui-datepicker-next-hover {
	right: 1px;
}
.ui-datepicker .ui-datepicker-prev span,
.ui-datepicker .ui-datepicker-next span {
	display: block;
	position: absolute;
	left: 50%;
	margin-left: -8px;
	top: 50%;
	margin-top: -8px;
}
.ui-datepicker .ui-datepicker-title {
	margin: 0 2.3em;
	line-height: 1.8em;
	text-align: center;
}
.ui-datepicker .ui-datepicker-title select {
	font-size: 1em;
	margin: 1px 0;
}
.ui-datepicker select.ui-datepicker-month,
.ui-datepicker select.ui-datepicker-year {
	width: 45%;
}
.ui-datepicker table {
	width: 100%;
	font-size: .9em;
	border-collapse: collapse;
	margin: 0 0 .4em;
}
.ui-datepicker th {
	padding: .7em .3em;
	text-align: center;
	font-weight: bold;
	border: 0;
}
.ui-datepicker td {
	border: 0;
	padding: 1px;
}
.ui-datepicker td span,
.ui-datepicker td a {
	display: block;
	padding: .2em;
	text-align: right;
	text-decoration: none;
}
.ui-datepicker .ui-datepicker-buttonpane {
	background-image: none;
	margin: .7em 0 0 0;
	padding: 0 .2em;
	border-left: 0;
	border-right: 0;
	border-bottom: 0;
}
.ui-datepicker .ui-datepicker-buttonpane button {
	float: right;
	margin: .5em .2em .4em;
	cursor: pointer;
	padding: .2em .6em .3em .6em;
	width: auto;
	overflow: visible;
}
.ui-datepicker .ui-datepicker-buttonpane button.ui-datepicker-current {
	float: left;
}


.ui-datepicker.ui-datepicker-multi {
	width: auto;
}
.ui-datepicker-multi .ui-datepicker-group {
	float: left;
}
.ui-datepicker-multi .ui-datepicker-group table {
	width: 95%;
	margin: 0 auto .4em;
}
.ui-datepicker-multi-2 .ui-datepicker-group {
	width: 50%;
}
.ui-datepicker-multi-3 .ui-datepicker-group {
	width: 33.3%;
}
.ui-datepicker-multi-4 .ui-datepicker-group {
	width: 25%;
}
.ui-datepicker-multi .ui-datepicker-group-last .ui-datepicker-header,
.ui-datepicker-multi .ui-datepicker-group-middle .ui-datepicker-header {
	border-left-width: 0;
}
.ui-datepicker-multi .ui-datepicker-buttonpane {
	clear: left;
}
.ui-datepicker-row-break {
	clear: both;
	width: 100%;
	font-size: 0;
}


.ui-datepicker-rtl {
	direction: rtl;
}
.ui-datepicker-rtl .ui-datepicker-prev {
	right: 2px;
	left: auto;
}
.ui-datepicker-rtl .ui-datepicker-next {
	left: 2px;
	right: auto;
}
.ui-datepicker-rtl .ui-datepicker-prev:hover {
	right: 1px;
	left: auto;
}
.ui-datepicker-rtl .ui-datepicker-next:hover {
	left: 1px;
	right: auto;
}
.ui-datepicker-rtl .ui-datepicker-buttonpane {
	clear: right;
}
.ui-datepicker-rtl .ui-datepicker-buttonpane button {
	float: left;
}
.ui-datepicker-rtl .ui-datepicker-buttonpane button.ui-datepicker-current,
.ui-datepicker-rtl .ui-datepicker-group {
	float: right;
}
.ui-datepicker-rtl .ui-datepicker-group-last .ui-datepicker-header,
.ui-datepicker-rtl .ui-datepicker-group-middle .ui-datepicker-header {
	border-right-width: 0;
	border-left-width: 1px;
}


.ui-datepicker .ui-icon {
	display: block;
	text-indent: -99999px;
	overflow: hidden;
	background-repeat: no-repeat;
	left: .5em;
	top: .3em;
}
.ui-dialog {
	position: absolute;
	top: 0;
	left: 0;
	padding: .2em;
	outline: 0;
}
.ui-dialog .ui-dialog-titlebar {
	padding: .4em 1em;
	position: relative;
}
.ui-dialog .ui-dialog-title {
	float: left;
	margin: .1em 0;
	white-space: nowrap;
	width: 90%;
	overflow: hidden;
	text-overflow: ellipsis;
}
.ui-dialog .ui-dialog-titlebar-close {
	position: absolute;
	right: .3em;
	top: 50%;
	width: 20px;
	margin: -10px 0 0 0;
	padding: 1px;
	height: 20px;
}
.ui-dialog .ui-dialog-content {
	position: relative;
	border: 0;
	padding: .5em 1em;
	background: none;
	overflow: auto;
}
.ui-dialog .ui-dialog-buttonpane {
	text-align: left;
	border-width: 1px 0 0 0;
	background-image: none;
	margin-top: .5em;
	padding: .3em 1em .5em .4em;
}
.ui-dialog .ui-dialog-buttonpane .ui-dialog-buttonset {
	float: right;
}
.ui-dialog .ui-dialog-buttonpane button {
	margin: .5em .4em .5em 0;
	cursor: pointer;
}
.ui-dialog .ui-resizable-n {
	height: 2px;
	top: 0;
}
.ui-dialog .ui-resizable-e {
	width: 2px;
	right: 0;
}
.ui-dialog .ui-resizable-s {
	height: 2px;
	bottom: 0;
}
.ui-dialog .ui-resizable-w {
	width: 2px;
	left: 0;
}
.ui-dialog .ui-resizable-se,
.ui-dialog .ui-resizable-sw,
.ui-dialog .ui-resizable-ne,
.ui-dialog .ui-resizable-nw {
	width: 7px;
	height: 7px;
}
.ui-dialog .ui-resizable-se {
	right: 0;
	bottom: 0;
}
.ui-dialog .ui-resizable-sw {
	left: 0;
	bottom: 0;
}
.ui-dialog .ui-resizable-ne {
	right: 0;
	top: 0;
}
.ui-dialog .ui-resizable-nw {
	left: 0;
	top: 0;
}
.ui-draggable .ui-dialog-titlebar {
	cursor: move;
}
.ui-draggable-handle {
	-ms-touch-action: none;
	touch-action: none;
}
.ui-resizable {
	position: relative;
}
.ui-resizable-handle {
	position: absolute;
	font-size: 0.1px;
	display: block;
	-ms-touch-action: none;
	touch-action: none;
}
.ui-resizable-disabled .ui-resizable-handle,
.ui-resizable-autohide .ui-resizable-handle {
	display: none;
}
.ui-resizable-n {
	cursor: n-resize;
	height: 7px;
	width: 100%;
	top: -5px;
	left: 0;
}
.ui-resizable-s {
	cursor: s-resize;
	height: 7px;
	width: 100%;
	bottom: -5px;
	left: 0;
}
.ui-resizable-e {
	cursor: e-resize;
	width: 7px;
	right: -5px;
	top: 0;
	height: 100%;
}
.ui-resizable-w {
	cursor: w-resize;
	width: 7px;
	left: -5px;
	top: 0;
	height: 100%;
}
.ui-resizable-se {
	cursor: se-resize;
	width: 12px;
	height: 12px;
	right: 1px;
	bottom: 1px;
}
.ui-resizable-sw {
	cursor: sw-resize;
	width: 9px;
	height: 9px;
	left: -5px;
	bottom: -5px;
}
.ui-resizable-nw {
	cursor: nw-resize;
	width: 9px;
	height: 9px;
	left: -5px;
	top: -5px;
}
.ui-resizable-ne {
	cursor: ne-resize;
	width: 9px;
	height: 9px;
	right: -5px;
	top: -5px;
}
.ui-progressbar {
	height: 2em;
	text-align: left;
	overflow: hidden;
}
.ui-progressbar .ui-progressbar-value {
	margin: -1px;
	height: 100%;
}
.ui-progressbar .ui-progressbar-overlay {
	background: url("data:image/gif;base64,R0lGODlhKAAoAIABAAAAAP///yH/C05FVFNDQVBFMi4wAwEAAAAh+QQJAQABACwAAAAAKAAoAAACkYwNqXrdC52DS06a7MFZI+4FHBCKoDeWKXqymPqGqxvJrXZbMx7Ttc+w9XgU2FB3lOyQRWET2IFGiU9m1frDVpxZZc6bfHwv4c1YXP6k1Vdy292Fb6UkuvFtXpvWSzA+HycXJHUXiGYIiMg2R6W459gnWGfHNdjIqDWVqemH2ekpObkpOlppWUqZiqr6edqqWQAAIfkECQEAAQAsAAAAACgAKAAAApSMgZnGfaqcg1E2uuzDmmHUBR8Qil95hiPKqWn3aqtLsS18y7G1SzNeowWBENtQd+T1JktP05nzPTdJZlR6vUxNWWjV+vUWhWNkWFwxl9VpZRedYcflIOLafaa28XdsH/ynlcc1uPVDZxQIR0K25+cICCmoqCe5mGhZOfeYSUh5yJcJyrkZWWpaR8doJ2o4NYq62lAAACH5BAkBAAEALAAAAAAoACgAAAKVDI4Yy22ZnINRNqosw0Bv7i1gyHUkFj7oSaWlu3ovC8GxNso5fluz3qLVhBVeT/Lz7ZTHyxL5dDalQWPVOsQWtRnuwXaFTj9jVVh8pma9JjZ4zYSj5ZOyma7uuolffh+IR5aW97cHuBUXKGKXlKjn+DiHWMcYJah4N0lYCMlJOXipGRr5qdgoSTrqWSq6WFl2ypoaUAAAIfkECQEAAQAsAAAAACgAKAAAApaEb6HLgd/iO7FNWtcFWe+ufODGjRfoiJ2akShbueb0wtI50zm02pbvwfWEMWBQ1zKGlLIhskiEPm9R6vRXxV4ZzWT2yHOGpWMyorblKlNp8HmHEb/lCXjcW7bmtXP8Xt229OVWR1fod2eWqNfHuMjXCPkIGNileOiImVmCOEmoSfn3yXlJWmoHGhqp6ilYuWYpmTqKUgAAIfkECQEAAQAsAAAAACgAKAAAApiEH6kb58biQ3FNWtMFWW3eNVcojuFGfqnZqSebuS06w5V80/X02pKe8zFwP6EFWOT1lDFk8rGERh1TTNOocQ61Hm4Xm2VexUHpzjymViHrFbiELsefVrn6XKfnt2Q9G/+Xdie499XHd2g4h7ioOGhXGJboGAnXSBnoBwKYyfioubZJ2Hn0RuRZaflZOil56Zp6iioKSXpUAAAh+QQJAQABACwAAAAAKAAoAAACkoQRqRvnxuI7kU1a1UU5bd5tnSeOZXhmn5lWK3qNTWvRdQxP8qvaC+/yaYQzXO7BMvaUEmJRd3TsiMAgswmNYrSgZdYrTX6tSHGZO73ezuAw2uxuQ+BbeZfMxsexY35+/Qe4J1inV0g4x3WHuMhIl2jXOKT2Q+VU5fgoSUI52VfZyfkJGkha6jmY+aaYdirq+lQAACH5BAkBAAEALAAAAAAoACgAAAKWBIKpYe0L3YNKToqswUlvznigd4wiR4KhZrKt9Upqip61i9E3vMvxRdHlbEFiEXfk9YARYxOZZD6VQ2pUunBmtRXo1Lf8hMVVcNl8JafV38aM2/Fu5V16Bn63r6xt97j09+MXSFi4BniGFae3hzbH9+hYBzkpuUh5aZmHuanZOZgIuvbGiNeomCnaxxap2upaCZsq+1kAACH5BAkBAAEALAAAAAAoACgAAAKXjI8By5zf4kOxTVrXNVlv1X0d8IGZGKLnNpYtm8Lr9cqVeuOSvfOW79D9aDHizNhDJidFZhNydEahOaDH6nomtJjp1tutKoNWkvA6JqfRVLHU/QUfau9l2x7G54d1fl995xcIGAdXqMfBNadoYrhH+Mg2KBlpVpbluCiXmMnZ2Sh4GBqJ+ckIOqqJ6LmKSllZmsoq6wpQAAAh+QQJAQABACwAAAAAKAAoAAAClYx/oLvoxuJDkU1a1YUZbJ59nSd2ZXhWqbRa2/gF8Gu2DY3iqs7yrq+xBYEkYvFSM8aSSObE+ZgRl1BHFZNr7pRCavZ5BW2142hY3AN/zWtsmf12p9XxxFl2lpLn1rseztfXZjdIWIf2s5dItwjYKBgo9yg5pHgzJXTEeGlZuenpyPmpGQoKOWkYmSpaSnqKileI2FAAACH5BAkBAAEALAAAAAAoACgAAAKVjB+gu+jG4kORTVrVhRlsnn2dJ3ZleFaptFrb+CXmO9OozeL5VfP99HvAWhpiUdcwkpBH3825AwYdU8xTqlLGhtCosArKMpvfa1mMRae9VvWZfeB2XfPkeLmm18lUcBj+p5dnN8jXZ3YIGEhYuOUn45aoCDkp16hl5IjYJvjWKcnoGQpqyPlpOhr3aElaqrq56Bq7VAAAOw==");
	height: 100%;
	filter: alpha(opacity=25); 
	opacity: 0.25;
}
.ui-progressbar-indeterminate .ui-progressbar-value {
	background-image: none;
}
.ui-selectable {
	-ms-touch-action: none;
	touch-action: none;
}
.ui-selectable-helper {
	position: absolute;
	z-index: 100;
	border: 1px dotted black;
}
.ui-selectmenu-menu {
	padding: 0;
	margin: 0;
	position: absolute;
	top: 0;
	left: 0;
	display: none;
}
.ui-selectmenu-menu .ui-menu {
	overflow: auto;
	overflow-x: hidden;
	padding-bottom: 1px;
}
.ui-selectmenu-menu .ui-menu .ui-selectmenu-optgroup {
	font-size: 1em;
	font-weight: bold;
	line-height: 1.5;
	padding: 2px 0.4em;
	margin: 0.5em 0 0 0;
	height: auto;
	border: 0;
}
.ui-selectmenu-open {
	display: block;
}
.ui-selectmenu-text {
	display: block;
	margin-right: 20px;
	overflow: hidden;
	text-overflow: ellipsis;
}
.ui-selectmenu-button.ui-button {
	text-align: left;
	white-space: nowrap;
	width: 14em;
}
.ui-selectmenu-icon.ui-icon {
	float: right;
	margin-top: 0;
}
.ui-slider {
	position: relative;
	text-align: left;
}
.ui-slider .ui-slider-handle {
	position: absolute;
	z-index: 2;
	width: 1.2em;
	height: 1.2em;
	cursor: default;
	-ms-touch-action: none;
	touch-action: none;
}
.ui-slider .ui-slider-range {
	position: absolute;
	z-index: 1;
	font-size: .7em;
	display: block;
	border: 0;
	background-position: 0 0;
}


.ui-slider.ui-state-disabled .ui-slider-handle,
.ui-slider.ui-state-disabled .ui-slider-range {
	filter: inherit;
}

.ui-slider-horizontal {
	height: .8em;
}
.ui-slider-horizontal .ui-slider-handle {
	top: -.3em;
	margin-left: -.6em;
}
.ui-slider-horizontal .ui-slider-range {
	top: 0;
	height: 100%;
}
.ui-slider-horizontal .ui-slider-range-min {
	left: 0;
}
.ui-slider-horizontal .ui-slider-range-max {
	right: 0;
}

.ui-slider-vertical {
	width: .8em;
	height: 100px;
}
.ui-slider-vertical .ui-slider-handle {
	left: -.3em;
	margin-left: 0;
	margin-bottom: -.6em;
}
.ui-slider-vertical .ui-slider-range {
	left: 0;
	width: 100%;
}
.ui-slider-vertical .ui-slider-range-min {
	bottom: 0;
}
.ui-slider-vertical .ui-slider-range-max {
	top: 0;
}
.ui-sortable-handle {
	-ms-touch-action: none;
	touch-action: none;
}
.ui-spinner {
	position: relative;
	display: inline-block;
	overflow: hidden;
	padding: 0;
	vertical-align: middle;
}
.ui-spinner-input {
	border: none;
	background: none;
	color: inherit;
	padding: .222em 0;
	margin: .2em 0;
	vertical-align: middle;
	margin-left: .4em;
	margin-right: 2em;
}
.ui-spinner-button {
	width: 1.6em;
	height: 50%;
	font-size: .5em;
	padding: 0;
	margin: 0;
	text-align: center;
	position: absolute;
	cursor: default;
	display: block;
	overflow: hidden;
	right: 0;
}

.ui-spinner a.ui-spinner-button {
	border-top-style: none;
	border-bottom-style: none;
	border-right-style: none;
}
.ui-spinner-up {
	top: 0;
}
.ui-spinner-down {
	bottom: 0;
}
.ui-tabs {
	position: relative;
	padding: .2em;
}
.ui-tabs .ui-tabs-nav {
	margin: 0;
	padding: .2em .2em 0;
}
.ui-tabs .ui-tabs-nav li {
	list-style: none;
	float: left;
	position: relative;
	top: 0;
	margin: 1px .2em 0 0;
	border-bottom-width: 0;
	padding: 0;
	white-space: nowrap;
}
.ui-tabs .ui-tabs-nav .ui-tabs-anchor {
	float: left;
	padding: .5em 1em;
	text-decoration: none;
}
.ui-tabs .ui-tabs-nav li.ui-tabs-active {
	margin-bottom: -1px;
	padding-bottom: 1px;
}
.ui-tabs .ui-tabs-nav li.ui-tabs-active .ui-tabs-anchor,
.ui-tabs .ui-tabs-nav li.ui-state-disabled .ui-tabs-anchor,
.ui-tabs .ui-tabs-nav li.ui-tabs-loading .ui-tabs-anchor {
	cursor: text;
}
.ui-tabs-collapsible .ui-tabs-nav li.ui-tabs-active .ui-tabs-anchor {
	cursor: pointer;
}
.ui-tabs .ui-tabs-panel {
	display: block;
	border-width: 0;
	padding: 1em 1.4em;
	background: none;
}
.ui-tooltip {
	padding: 8px;
	position: absolute;
	z-index: 9999;
	max-width: 300px;
}
body .ui-tooltip {
	border-width: 2px;
}

.ui-widget {
	font-family: Arial,Helvetica,sans-serif;
	font-size: 1em;
}
.ui-widget .ui-widget {
	font-size: 1em;
}
.ui-widget input,
.ui-widget select,
.ui-widget textarea,
.ui-widget button {
	font-family: Arial,Helvetica,sans-serif;
	font-size: 1em;
}
.ui-widget.ui-widget-content {
	border: 1px solid #c5c5c5;
}
.ui-widget-content {
	border: 1px solid #dddddd;
	background: #ffffff;
	color: #333333;
}
.ui-widget-content a {
	color: #333333;
}
.ui-widget-header {
	border: 1px solid #dddddd;
	background: #e9e9e9;
	color: #333333;
	font-weight: bold;
}
.ui-widget-header a {
	color: #333333;
}


.ui-state-default,
.ui-widget-content .ui-state-default,
.ui-widget-header .ui-state-default,
.ui-button,


html .ui-button.ui-state-disabled:hover,
html .ui-button.ui-state-disabled:active {
	border: 1px solid #c5c5c5;
	background: #f6f6f6;
	font-weight: normal;
	color: #454545;
}
.ui-state-default a,
.ui-state-default a:link,
.ui-state-default a:visited,
a.ui-button,
a:link.ui-button,
a:visited.ui-button,
.ui-button {
	color: #454545;
	text-decoration: none;
}
.ui-state-hover,
.ui-widget-content .ui-state-hover,
.ui-widget-header .ui-state-hover,
.ui-state-focus,
.ui-widget-content .ui-state-focus,
.ui-widget-header .ui-state-focus,
.ui-button:hover,
.ui-button:focus {
	border: 1px solid #cccccc;
	background: #ededed;
	font-weight: normal;
	color: #2b2b2b;
}
.ui-state-hover a,
.ui-state-hover a:hover,
.ui-state-hover a:link,
.ui-state-hover a:visited,
.ui-state-focus a,
.ui-state-focus a:hover,
.ui-state-focus a:link,
.ui-state-focus a:visited,
a.ui-button:hover,
a.ui-button:focus {
	color: #2b2b2b;
	text-decoration: none;
}

.ui-visual-focus {
	box-shadow: 0 0 3px 1px rgb(94, 158, 214);
}
.ui-state-active,
.ui-widget-content .ui-state-active,
.ui-widget-header .ui-state-active,
a.ui-button:active,
.ui-button:active,
.ui-button.ui-state-active:hover {
	border: 1px solid #003eff;
	background: #007fff;
	font-weight: normal;
	color: #ffffff;
}
.ui-icon-background,
.ui-state-active .ui-icon-background {
	border: #003eff;
	background-color: #ffffff;
}
.ui-state-active a,
.ui-state-active a:link,
.ui-state-active a:visited {
	color: #ffffff;
	text-decoration: none;
}


.ui-state-highlight,
.ui-widget-content .ui-state-highlight,
.ui-widget-header .ui-state-highlight {
	border: 1px solid #dad55e;
	background: #fffa90;
	color: #777620;
}
.ui-state-checked {
	border: 1px solid #dad55e;
	background: #fffa90;
}
.ui-state-highlight a,
.ui-widget-content .ui-state-highlight a,
.ui-widget-header .ui-state-highlight a {
	color: #777620;
}
.ui-state-error,
.ui-widget-content .ui-state-error,
.ui-widget-header .ui-state-error {
	border: 1px solid #f1a899;
	background: #fddfdf;
	color: #5f3f3f;
}
.ui-state-error a,
.ui-widget-content .ui-state-error a,
.ui-widget-header .ui-state-error a {
	color: #5f3f3f;
}
.ui-state-error-text,
.ui-widget-content .ui-state-error-text,
.ui-widget-header .ui-state-error-text {
	color: #5f3f3f;
}
.ui-priority-primary,
.ui-widget-content .ui-priority-primary,
.ui-widget-header .ui-priority-primary {
	font-weight: bold;
}
.ui-priority-secondary,
.ui-widget-content .ui-priority-secondary,
.ui-widget-header .ui-priority-secondary {
	opacity: .7;
	filter:Alpha(Opacity=70); 
	font-weight: normal;
}
.ui-state-disabled,
.ui-widget-content .ui-state-disabled,
.ui-widget-header .ui-state-disabled {
	opacity: .35;
	filter:Alpha(Opacity=35); 
	background-image: none;
}
.ui-state-disabled .ui-icon {
	filter:Alpha(Opacity=35); 
}




.ui-icon {
	width: 16px;
	height: 16px;
}
.ui-icon,
.ui-widget-content .ui-icon {
	background-image: url("images/ui-icons_444444_256x240.png");
}
.ui-widget-header .ui-icon {
	background-image: url("images/ui-icons_444444_256x240.png");
}
.ui-state-hover .ui-icon,
.ui-state-focus .ui-icon,
.ui-button:hover .ui-icon,
.ui-button:focus .ui-icon {
	background-image: url("images/ui-icons_555555_256x240.png");
}
.ui-state-active .ui-icon,
.ui-button:active .ui-icon {
	background-image: url("images/ui-icons_ffffff_256x240.png");
}
.ui-state-highlight .ui-icon,
.ui-button .ui-state-highlight.ui-icon {
	background-image: url("images/ui-icons_777620_256x240.png");
}
.ui-state-error .ui-icon,
.ui-state-error-text .ui-icon {
	background-image: url("images/ui-icons_cc0000_256x240.png");
}
.ui-button .ui-icon {
	background-image: url("images/ui-icons_777777_256x240.png");
}


.ui-icon-blank { background-position: 16px 16px; }
.ui-icon-caret-1-n { background-position: 0 0; }
.ui-icon-caret-1-ne { background-position: -16px 0; }
.ui-icon-caret-1-e { background-position: -32px 0; }
.ui-icon-caret-1-se { background-position: -48px 0; }
.ui-icon-caret-1-s { background-position: -65px 0; }
.ui-icon-caret-1-sw { background-position: -80px 0; }
.ui-icon-caret-1-w { background-position: -96px 0; }
.ui-icon-caret-1-nw { background-position: -112px 0; }
.ui-icon-caret-2-n-s { background-position: -128px 0; }
.ui-icon-caret-2-e-w { background-position: -144px 0; }
.ui-icon-triangle-1-n { background-position: 0 -16px; }
.ui-icon-triangle-1-ne { background-position: -16px -16px; }
.ui-icon-triangle-1-e { background-position: -32px -16px; }
.ui-icon-triangle-1-se { background-position: -48px -16px; }
.ui-icon-triangle-1-s { background-position: -65px -16px; }
.ui-icon-triangle-1-sw { background-position: -80px -16px; }
.ui-icon-triangle-1-w { background-position: -96px -16px; }
.ui-icon-triangle-1-nw { background-position: -112px -16px; }
.ui-icon-triangle-2-n-s { background-position: -128px -16px; }
.ui-icon-triangle-2-e-w { background-position: -144px -16px; }
.ui-icon-arrow-1-n { background-position: 0 -32px; }
.ui-icon-arrow-1-ne { background-position: -16px -32px; }
.ui-icon-arrow-1-e { background-position: -32px -32px; }
.ui-icon-arrow-1-se { background-position: -48px -32px; }
.ui-icon-arrow-1-s { background-position: -65px -32px; }
.ui-icon-arrow-1-sw { background-position: -80px -32px; }
.ui-icon-arrow-1-w { background-position: -96px -32px; }
.ui-icon-arrow-1-nw { background-position: -112px -32px; }
.ui-icon-arrow-2-n-s { background-position: -128px -32px; }
.ui-icon-arrow-2-ne-sw { background-position: -144px -32px; }
.ui-icon-arrow-2-e-w { background-position: -160px -32px; }
.ui-icon-arrow-2-se-nw { background-position: -176px -32px; }
.ui-icon-arrowstop-1-n { background-position: -192px -32px; }
.ui-icon-arrowstop-1-e { background-position: -208px -32px; }
.ui-icon-arrowstop-1-s { background-position: -224px -32px; }
.ui-icon-arrowstop-1-w { background-position: -240px -32px; }
.ui-icon-arrowthick-1-n { background-position: 1px -48px; }
.ui-icon-arrowthick-1-ne { background-position: -16px -48px; }
.ui-icon-arrowthick-1-e { background-position: -32px -48px; }
.ui-icon-arrowthick-1-se { background-position: -48px -48px; }
.ui-icon-arrowthick-1-s { background-position: -64px -48px; }
.ui-icon-arrowthick-1-sw { background-position: -80px -48px; }
.ui-icon-arrowthick-1-w { background-position: -96px -48px; }
.ui-icon-arrowthick-1-nw { background-position: -112px -48px; }
.ui-icon-arrowthick-2-n-s { background-position: -128px -48px; }
.ui-icon-arrowthick-2-ne-sw { background-position: -144px -48px; }
.ui-icon-arrowthick-2-e-w { background-position: -160px -48px; }
.ui-icon-arrowthick-2-se-nw { background-position: -176px -48px; }
.ui-icon-arrowthickstop-1-n { background-position: -192px -48px; }
.ui-icon-arrowthickstop-1-e { background-position: -208px -48px; }
.ui-icon-arrowthickstop-1-s { background-position: -224px -48px; }
.ui-icon-arrowthickstop-1-w { background-position: -240px -48px; }
.ui-icon-arrowreturnthick-1-w { background-position: 0 -64px; }
.ui-icon-arrowreturnthick-1-n { background-position: -16px -64px; }
.ui-icon-arrowreturnthick-1-e { background-position: -32px -64px; }
.ui-icon-arrowreturnthick-1-s { background-position: -48px -64px; }
.ui-icon-arrowreturn-1-w { background-position: -64px -64px; }
.ui-icon-arrowreturn-1-n { background-position: -80px -64px; }
.ui-icon-arrowreturn-1-e { background-position: -96px -64px; }
.ui-icon-arrowreturn-1-s { background-position: -112px -64px; }
.ui-icon-arrowrefresh-1-w { background-position: -128px -64px; }
.ui-icon-arrowrefresh-1-n { background-position: -144px -64px; }
.ui-icon-arrowrefresh-1-e { background-position: -160px -64px; }
.ui-icon-arrowrefresh-1-s { background-position: -176px -64px; }
.ui-icon-arrow-4 { background-position: 0 -80px; }
.ui-icon-arrow-4-diag { background-position: -16px -80px; }
.ui-icon-extlink { background-position: -32px -80px; }
.ui-icon-newwin { background-position: -48px -80px; }
.ui-icon-refresh { background-position: -64px -80px; }
.ui-icon-shuffle { background-position: -80px -80px; }
.ui-icon-transfer-e-w { background-position: -96px -80px; }
.ui-icon-transferthick-e-w { background-position: -112px -80px; }
.ui-icon-folder-collapsed { background-position: 0 -96px; }
.ui-icon-folder-open { background-position: -16px -96px; }
.ui-icon-document { background-position: -32px -96px; }
.ui-icon-document-b { background-position: -48px -96px; }
.ui-icon-note { background-position: -64px -96px; }
.ui-icon-mail-closed { background-position: -80px -96px; }
.ui-icon-mail-open { background-position: -96px -96px; }
.ui-icon-suitcase { background-position: -112px -96px; }
.ui-icon-comment { background-position: -128px -96px; }
.ui-icon-person { background-position: -144px -96px; }
.ui-icon-print { background-position: -160px -96px; }
.ui-icon-trash { background-position: -176px -96px; }
.ui-icon-locked { background-position: -192px -96px; }
.ui-icon-unlocked { background-position: -208px -96px; }
.ui-icon-bookmark { background-position: -224px -96px; }
.ui-icon-tag { background-position: -240px -96px; }
.ui-icon-home { background-position: 0 -112px; }
.ui-icon-flag { background-position: -16px -112px; }
.ui-icon-calendar { background-position: -32px -112px; }
.ui-icon-cart { background-position: -48px -112px; }
.ui-icon-pencil { background-position: -64px -112px; }
.ui-icon-clock { background-position: -80px -112px; }
.ui-icon-disk { background-position: -96px -112px; }
.ui-icon-calculator { background-position: -112px -112px; }
.ui-icon-zoomin { background-position: -128px -112px; }
.ui-icon-zoomout { background-position: -144px -112px; }
.ui-icon-search { background-position: -160px -112px; }
.ui-icon-wrench { background-position: -176px -112px; }
.ui-icon-gear { background-position: -192px -112px; }
.ui-icon-heart { background-position: -208px -112px; }
.ui-icon-star { background-position: -224px -112px; }
.ui-icon-link { background-position: -240px -112px; }
.ui-icon-cancel { background-position: 0 -128px; }
.ui-icon-plus { background-position: -16px -128px; }
.ui-icon-plusthick { background-position: -32px -128px; }
.ui-icon-minus { background-position: -48px -128px; }
.ui-icon-minusthick { background-position: -64px -128px; }
.ui-icon-close { background-position: -80px -128px; }
.ui-icon-closethick { background-position: -96px -128px; }
.ui-icon-key { background-position: -112px -128px; }
.ui-icon-lightbulb { background-position: -128px -128px; }
.ui-icon-scissors { background-position: -144px -128px; }
.ui-icon-clipboard { background-position: -160px -128px; }
.ui-icon-copy { background-position: -176px -128px; }
.ui-icon-contact { background-position: -192px -128px; }
.ui-icon-image { background-position: -208px -128px; }
.ui-icon-video { background-position: -224px -128px; }
.ui-icon-script { background-position: -240px -128px; }
.ui-icon-alert { background-position: 0 -144px; }
.ui-icon-info { background-position: -16px -144px; }
.ui-icon-notice { background-position: -32px -144px; }
.ui-icon-help { background-position: -48px -144px; }
.ui-icon-check { background-position: -64px -144px; }
.ui-icon-bullet { background-position: -80px -144px; }
.ui-icon-radio-on { background-position: -96px -144px; }
.ui-icon-radio-off { background-position: -112px -144px; }
.ui-icon-pin-w { background-position: -128px -144px; }
.ui-icon-pin-s { background-position: -144px -144px; }
.ui-icon-play { background-position: 0 -160px; }
.ui-icon-pause { background-position: -16px -160px; }
.ui-icon-seek-next { background-position: -32px -160px; }
.ui-icon-seek-prev { background-position: -48px -160px; }
.ui-icon-seek-end { background-position: -64px -160px; }
.ui-icon-seek-start { background-position: -80px -160px; }

.ui-icon-seek-first { background-position: -80px -160px; }
.ui-icon-stop { background-position: -96px -160px; }
.ui-icon-eject { background-position: -112px -160px; }
.ui-icon-volume-off { background-position: -128px -160px; }
.ui-icon-volume-on { background-position: -144px -160px; }
.ui-icon-power { background-position: 0 -176px; }
.ui-icon-signal-diag { background-position: -16px -176px; }
.ui-icon-signal { background-position: -32px -176px; }
.ui-icon-battery-0 { background-position: -48px -176px; }
.ui-icon-battery-1 { background-position: -64px -176px; }
.ui-icon-battery-2 { background-position: -80px -176px; }
.ui-icon-battery-3 { background-position: -96px -176px; }
.ui-icon-circle-plus { background-position: 0 -192px; }
.ui-icon-circle-minus { background-position: -16px -192px; }
.ui-icon-circle-close { background-position: -32px -192px; }
.ui-icon-circle-triangle-e { background-position: -48px -192px; }
.ui-icon-circle-triangle-s { background-position: -64px -192px; }
.ui-icon-circle-triangle-w { background-position: -80px -192px; }
.ui-icon-circle-triangle-n { background-position: -96px -192px; }
.ui-icon-circle-arrow-e { background-position: -112px -192px; }
.ui-icon-circle-arrow-s { background-position: -128px -192px; }
.ui-icon-circle-arrow-w { background-position: -144px -192px; }
.ui-icon-circle-arrow-n { background-position: -160px -192px; }
.ui-icon-circle-zoomin { background-position: -176px -192px; }
.ui-icon-circle-zoomout { background-position: -192px -192px; }
.ui-icon-circle-check { background-position: -208px -192px; }
.ui-icon-circlesmall-plus { background-position: 0 -208px; }
.ui-icon-circlesmall-minus { background-position: -16px -208px; }
.ui-icon-circlesmall-close { background-position: -32px -208px; }
.ui-icon-squaresmall-plus { background-position: -48px -208px; }
.ui-icon-squaresmall-minus { background-position: -64px -208px; }
.ui-icon-squaresmall-close { background-position: -80px -208px; }
.ui-icon-grip-dotted-vertical { background-position: 0 -224px; }
.ui-icon-grip-dotted-horizontal { background-position: -16px -224px; }
.ui-icon-grip-solid-vertical { background-position: -32px -224px; }
.ui-icon-grip-solid-horizontal { background-position: -48px -224px; }
.ui-icon-gripsmall-diagonal-se { background-position: -64px -224px; }
.ui-icon-grip-diagonal-se { background-position: -80px -224px; }





.ui-corner-all,
.ui-corner-top,
.ui-corner-left,
.ui-corner-tl {
	border-top-left-radius: 3px;
}
.ui-corner-all,
.ui-corner-top,
.ui-corner-right,
.ui-corner-tr {
	border-top-right-radius: 3px;
}
.ui-corner-all,
.ui-corner-bottom,
.ui-corner-left,
.ui-corner-bl {
	border-bottom-left-radius: 3px;
}
.ui-corner-all,
.ui-corner-bottom,
.ui-corner-right,
.ui-corner-br {
	border-bottom-right-radius: 3px;
}


.ui-widget-overlay {
	background: #aaaaaa;
	opacity: .3;
	filter: Alpha(Opacity=30); 
}
.ui-widget-shadow {
	-webkit-box-shadow: 0px 0px 5px #666666;
	box-shadow: 0px 0px 5px #666666;
} 
</style>
<script src='js/jquery-1.12.4.js'></script>
<script src='js/jquery-ui.js'></script>

<script>
$( function() {
	$( '#tabs-collect-classes-from-hierarchical-partition' ).tabs();
} );
</script>
<div id='tabs-collect-classes-from-hierarchical-partition'>
<ul>
<li><a href='#tab-collect-classes-from-hierarchical-partition-1'>n_signatures ≥ 134</a></li>
<li><a href='#tab-collect-classes-from-hierarchical-partition-2'>n_signatures ≥ 650</a></li>
<li><a href='#tab-collect-classes-from-hierarchical-partition-3'>n_signatures ≥ 2054</a></li>
</ul>
<div id='tab-collect-classes-from-hierarchical-partition-1'>
<pre><code class="r">collect_classes(res_rh, merge_node = merge_node_param(min_n_signatures = 134))
</code></pre>

<p><img src="figure_cola/tab-collect-classes-from-hierarchical-partition-1-1.png" alt="plot of chunk tab-collect-classes-from-hierarchical-partition-1"/></p>

</div>
<div id='tab-collect-classes-from-hierarchical-partition-2'>
<pre><code class="r">collect_classes(res_rh, merge_node = merge_node_param(min_n_signatures = 650))
</code></pre>

<p><img src="figure_cola/tab-collect-classes-from-hierarchical-partition-2-1.png" alt="plot of chunk tab-collect-classes-from-hierarchical-partition-2"/></p>

</div>
<div id='tab-collect-classes-from-hierarchical-partition-3'>
<pre><code class="r">collect_classes(res_rh, merge_node = merge_node_param(min_n_signatures = 2054))
</code></pre>

<p><img src="figure_cola/tab-collect-classes-from-hierarchical-partition-3-1.png" alt="plot of chunk tab-collect-classes-from-hierarchical-partition-3"/></p>

</div>
</div>

Following shows the table of the partitions (You need to click the **show/hide
code output** link to see it).


<script>
$( function() {
	$( '#tabs-get-classes-from-hierarchical-partition' ).tabs();
} );
</script>
<div id='tabs-get-classes-from-hierarchical-partition'>
<ul>
<li><a href='#tab-get-classes-from-hierarchical-partition-1'>n_signatures ≥ 134</a></li>
<li><a href='#tab-get-classes-from-hierarchical-partition-2'>n_signatures ≥ 650</a></li>
<li><a href='#tab-get-classes-from-hierarchical-partition-3'>n_signatures ≥ 2054</a></li>
</ul>

<div id='tab-get-classes-from-hierarchical-partition-1'>
<p><a id='tab-get-classes-from-hierarchical-partition-1-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">get_classes(res_rh, merge_node = merge_node_param(min_n_signatures = 134))
</code></pre>

<pre><code>#&gt; sample_39 sample_40 sample_42 sample_47 sample_48 sample_49 sample_41 sample_43 sample_44 sample_45 
#&gt;      &quot;03&quot;      &quot;03&quot;     &quot;013&quot;     &quot;012&quot;     &quot;011&quot;      &quot;03&quot;     &quot;011&quot;     &quot;012&quot;     &quot;012&quot;     &quot;012&quot; 
#&gt; sample_46 sample_70 sample_71 sample_72 sample_68 sample_69 sample_67 sample_55 sample_56 sample_59 
#&gt;     &quot;012&quot;     &quot;012&quot;     &quot;013&quot;     &quot;013&quot;     &quot;011&quot;     &quot;011&quot;     &quot;013&quot;      &quot;03&quot;      &quot;03&quot;     &quot;012&quot; 
#&gt; sample_52 sample_53 sample_51 sample_50 sample_54 sample_57 sample_58 sample_60 sample_61 sample_65 
#&gt;     &quot;022&quot;     &quot;021&quot;     &quot;021&quot;     &quot;021&quot;     &quot;012&quot;     &quot;022&quot;     &quot;022&quot;     &quot;012&quot;     &quot;022&quot;     &quot;022&quot; 
#&gt; sample_66 sample_63 sample_64 sample_62  sample_1  sample_2  sample_3  sample_4  sample_5  sample_6 
#&gt;     &quot;013&quot;     &quot;022&quot;     &quot;022&quot;     &quot;022&quot;      &quot;03&quot;     &quot;013&quot;      &quot;03&quot;      &quot;03&quot;     &quot;011&quot;      &quot;03&quot; 
#&gt;  sample_7  sample_8  sample_9 sample_10 sample_11 sample_12 sample_13 sample_14 sample_15 sample_16 
#&gt;      &quot;03&quot;      &quot;03&quot;     &quot;013&quot;     &quot;013&quot;     &quot;013&quot;     &quot;021&quot;     &quot;011&quot;     &quot;013&quot;     &quot;011&quot;     &quot;012&quot; 
#&gt; sample_17 sample_18 sample_19 sample_20 sample_21 sample_22 sample_23 sample_24 sample_25 sample_26 
#&gt;     &quot;011&quot;     &quot;012&quot;     &quot;012&quot;     &quot;011&quot;     &quot;011&quot;     &quot;021&quot;      &quot;03&quot;     &quot;011&quot;     &quot;021&quot;     &quot;012&quot; 
#&gt; sample_27 sample_34 sample_35 sample_36 sample_37 sample_38 sample_28 sample_29 sample_30 sample_31 
#&gt;      &quot;03&quot;     &quot;022&quot;     &quot;022&quot;     &quot;021&quot;     &quot;021&quot;     &quot;021&quot;     &quot;021&quot;     &quot;013&quot;     &quot;021&quot;     &quot;022&quot; 
#&gt; sample_32 sample_33 
#&gt;     &quot;021&quot;     &quot;021&quot;
</code></pre>

<script>
$('#tab-get-classes-from-hierarchical-partition-1-a').parent().next().next().hide();
$('#tab-get-classes-from-hierarchical-partition-1-a').click(function(){
  $('#tab-get-classes-from-hierarchical-partition-1-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-get-classes-from-hierarchical-partition-2'>
<p><a id='tab-get-classes-from-hierarchical-partition-2-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">get_classes(res_rh, merge_node = merge_node_param(min_n_signatures = 650))
</code></pre>

<pre><code>#&gt; sample_39 sample_40 sample_42 sample_47 sample_48 sample_49 sample_41 sample_43 sample_44 sample_45 
#&gt;      &quot;03&quot;      &quot;03&quot;     &quot;013&quot;     &quot;012&quot;     &quot;011&quot;      &quot;03&quot;     &quot;011&quot;     &quot;012&quot;     &quot;012&quot;     &quot;012&quot; 
#&gt; sample_46 sample_70 sample_71 sample_72 sample_68 sample_69 sample_67 sample_55 sample_56 sample_59 
#&gt;     &quot;012&quot;     &quot;012&quot;     &quot;013&quot;     &quot;013&quot;     &quot;011&quot;     &quot;011&quot;     &quot;013&quot;      &quot;03&quot;      &quot;03&quot;     &quot;012&quot; 
#&gt; sample_52 sample_53 sample_51 sample_50 sample_54 sample_57 sample_58 sample_60 sample_61 sample_65 
#&gt;      &quot;02&quot;      &quot;02&quot;      &quot;02&quot;      &quot;02&quot;     &quot;012&quot;      &quot;02&quot;      &quot;02&quot;     &quot;012&quot;      &quot;02&quot;      &quot;02&quot; 
#&gt; sample_66 sample_63 sample_64 sample_62  sample_1  sample_2  sample_3  sample_4  sample_5  sample_6 
#&gt;     &quot;013&quot;      &quot;02&quot;      &quot;02&quot;      &quot;02&quot;      &quot;03&quot;     &quot;013&quot;      &quot;03&quot;      &quot;03&quot;     &quot;011&quot;      &quot;03&quot; 
#&gt;  sample_7  sample_8  sample_9 sample_10 sample_11 sample_12 sample_13 sample_14 sample_15 sample_16 
#&gt;      &quot;03&quot;      &quot;03&quot;     &quot;013&quot;     &quot;013&quot;     &quot;013&quot;      &quot;02&quot;     &quot;011&quot;     &quot;013&quot;     &quot;011&quot;     &quot;012&quot; 
#&gt; sample_17 sample_18 sample_19 sample_20 sample_21 sample_22 sample_23 sample_24 sample_25 sample_26 
#&gt;     &quot;011&quot;     &quot;012&quot;     &quot;012&quot;     &quot;011&quot;     &quot;011&quot;      &quot;02&quot;      &quot;03&quot;     &quot;011&quot;      &quot;02&quot;     &quot;012&quot; 
#&gt; sample_27 sample_34 sample_35 sample_36 sample_37 sample_38 sample_28 sample_29 sample_30 sample_31 
#&gt;      &quot;03&quot;      &quot;02&quot;      &quot;02&quot;      &quot;02&quot;      &quot;02&quot;      &quot;02&quot;      &quot;02&quot;     &quot;013&quot;      &quot;02&quot;      &quot;02&quot; 
#&gt; sample_32 sample_33 
#&gt;      &quot;02&quot;      &quot;02&quot;
</code></pre>

<script>
$('#tab-get-classes-from-hierarchical-partition-2-a').parent().next().next().hide();
$('#tab-get-classes-from-hierarchical-partition-2-a').click(function(){
  $('#tab-get-classes-from-hierarchical-partition-2-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-get-classes-from-hierarchical-partition-3'>
<p><a id='tab-get-classes-from-hierarchical-partition-3-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">get_classes(res_rh, merge_node = merge_node_param(min_n_signatures = 2054))
</code></pre>

<pre><code>#&gt; sample_39 sample_40 sample_42 sample_47 sample_48 sample_49 sample_41 sample_43 sample_44 sample_45 
#&gt;      &quot;03&quot;      &quot;03&quot;      &quot;01&quot;      &quot;01&quot;      &quot;01&quot;      &quot;03&quot;      &quot;01&quot;      &quot;01&quot;      &quot;01&quot;      &quot;01&quot; 
#&gt; sample_46 sample_70 sample_71 sample_72 sample_68 sample_69 sample_67 sample_55 sample_56 sample_59 
#&gt;      &quot;01&quot;      &quot;01&quot;      &quot;01&quot;      &quot;01&quot;      &quot;01&quot;      &quot;01&quot;      &quot;01&quot;      &quot;03&quot;      &quot;03&quot;      &quot;01&quot; 
#&gt; sample_52 sample_53 sample_51 sample_50 sample_54 sample_57 sample_58 sample_60 sample_61 sample_65 
#&gt;      &quot;02&quot;      &quot;02&quot;      &quot;02&quot;      &quot;02&quot;      &quot;01&quot;      &quot;02&quot;      &quot;02&quot;      &quot;01&quot;      &quot;02&quot;      &quot;02&quot; 
#&gt; sample_66 sample_63 sample_64 sample_62  sample_1  sample_2  sample_3  sample_4  sample_5  sample_6 
#&gt;      &quot;01&quot;      &quot;02&quot;      &quot;02&quot;      &quot;02&quot;      &quot;03&quot;      &quot;01&quot;      &quot;03&quot;      &quot;03&quot;      &quot;01&quot;      &quot;03&quot; 
#&gt;  sample_7  sample_8  sample_9 sample_10 sample_11 sample_12 sample_13 sample_14 sample_15 sample_16 
#&gt;      &quot;03&quot;      &quot;03&quot;      &quot;01&quot;      &quot;01&quot;      &quot;01&quot;      &quot;02&quot;      &quot;01&quot;      &quot;01&quot;      &quot;01&quot;      &quot;01&quot; 
#&gt; sample_17 sample_18 sample_19 sample_20 sample_21 sample_22 sample_23 sample_24 sample_25 sample_26 
#&gt;      &quot;01&quot;      &quot;01&quot;      &quot;01&quot;      &quot;01&quot;      &quot;01&quot;      &quot;02&quot;      &quot;03&quot;      &quot;01&quot;      &quot;02&quot;      &quot;01&quot; 
#&gt; sample_27 sample_34 sample_35 sample_36 sample_37 sample_38 sample_28 sample_29 sample_30 sample_31 
#&gt;      &quot;03&quot;      &quot;02&quot;      &quot;02&quot;      &quot;02&quot;      &quot;02&quot;      &quot;02&quot;      &quot;02&quot;      &quot;01&quot;      &quot;02&quot;      &quot;02&quot; 
#&gt; sample_32 sample_33 
#&gt;      &quot;02&quot;      &quot;02&quot;
</code></pre>

<script>
$('#tab-get-classes-from-hierarchical-partition-3-a').parent().next().next().hide();
$('#tab-get-classes-from-hierarchical-partition-3-a').click(function(){
  $('#tab-get-classes-from-hierarchical-partition-3-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>
</div>



### Top rows heatmap

Heatmaps of the top rows:





```r
top_rows_heatmap(res_rh)
```

![plot of chunk top-rows-heatmap](figure_cola/top-rows-heatmap-1.png)

Top rows on each node:


```r
top_rows_overlap(res_rh, method = "upset")
```

![plot of chunk top-rows-overlap](figure_cola/top-rows-overlap-1.png)


### UMAP plot

UMAP plot which shows how samples are separated.




<script>
$( function() {
	$( '#tabs-dimension-reduction-by-depth' ).tabs();
} );
</script>
<div id='tabs-dimension-reduction-by-depth'>
<ul>
<li><a href='#tab-dimension-reduction-by-depth-1'>n_signatures ≥ 134</a></li>
<li><a href='#tab-dimension-reduction-by-depth-2'>n_signatures ≥ 650</a></li>
<li><a href='#tab-dimension-reduction-by-depth-3'>n_signatures ≥ 2054</a></li>
</ul>
<div id='tab-dimension-reduction-by-depth-1'>
<pre><code class="r">par(mfrow = c(1, 2))
dimension_reduction(res_rh, merge_node = merge_node_param(min_n_signatures = 134),
    method = &quot;UMAP&quot;, top_value_method = &quot;SD&quot;, top_n = 500, scale_rows = FALSE)
dimension_reduction(res_rh, merge_node = merge_node_param(min_n_signatures = 134),
    method = &quot;UMAP&quot;, top_value_method = &quot;ATC&quot;, top_n = 500, scale_rows = TRUE)
</code></pre>

<p><img src="figure_cola/tab-dimension-reduction-by-depth-1-1.png" title="plot of chunk tab-dimension-reduction-by-depth-1" alt="plot of chunk tab-dimension-reduction-by-depth-1" width="100%" /></p>

</div>
<div id='tab-dimension-reduction-by-depth-2'>
<pre><code class="r">par(mfrow = c(1, 2))
dimension_reduction(res_rh, merge_node = merge_node_param(min_n_signatures = 650),
    method = &quot;UMAP&quot;, top_value_method = &quot;SD&quot;, top_n = 500, scale_rows = FALSE)
dimension_reduction(res_rh, merge_node = merge_node_param(min_n_signatures = 650),
    method = &quot;UMAP&quot;, top_value_method = &quot;ATC&quot;, top_n = 500, scale_rows = TRUE)
</code></pre>

<p><img src="figure_cola/tab-dimension-reduction-by-depth-2-1.png" title="plot of chunk tab-dimension-reduction-by-depth-2" alt="plot of chunk tab-dimension-reduction-by-depth-2" width="100%" /></p>

</div>
<div id='tab-dimension-reduction-by-depth-3'>
<pre><code class="r">par(mfrow = c(1, 2))
dimension_reduction(res_rh, merge_node = merge_node_param(min_n_signatures = 2054),
    method = &quot;UMAP&quot;, top_value_method = &quot;SD&quot;, top_n = 500, scale_rows = FALSE)
dimension_reduction(res_rh, merge_node = merge_node_param(min_n_signatures = 2054),
    method = &quot;UMAP&quot;, top_value_method = &quot;ATC&quot;, top_n = 500, scale_rows = TRUE)
</code></pre>

<p><img src="figure_cola/tab-dimension-reduction-by-depth-3-1.png" title="plot of chunk tab-dimension-reduction-by-depth-3" alt="plot of chunk tab-dimension-reduction-by-depth-3" width="100%" /></p>

</div>
</div>




### Signature heatmap

Signatures on the heatmap are the union of all signatures found on every node
on the hierarchy. The number of k-means on rows are automatically selected by the function.




<script>
$( function() {
	$( '#tabs-get-signatures-from-hierarchical-partition' ).tabs();
} );
</script>
<div id='tabs-get-signatures-from-hierarchical-partition'>
<ul>
<li><a href='#tab-get-signatures-from-hierarchical-partition-1'>n_signatures ≥ 134</a></li>
<li><a href='#tab-get-signatures-from-hierarchical-partition-2'>n_signatures ≥ 650</a></li>
<li><a href='#tab-get-signatures-from-hierarchical-partition-3'>n_signatures ≥ 2054</a></li>
</ul>
<div id='tab-get-signatures-from-hierarchical-partition-1'>
<pre><code class="r">get_signatures(res_rh, merge_node = merge_node_param(min_n_signatures = 134))
</code></pre>

<p><img src="figure_cola/tab-get-signatures-from-hierarchical-partition-1-1.png" alt="plot of chunk tab-get-signatures-from-hierarchical-partition-1"/></p>

</div>
<div id='tab-get-signatures-from-hierarchical-partition-2'>
<pre><code class="r">get_signatures(res_rh, merge_node = merge_node_param(min_n_signatures = 650))
</code></pre>

<p><img src="figure_cola/tab-get-signatures-from-hierarchical-partition-2-1.png" alt="plot of chunk tab-get-signatures-from-hierarchical-partition-2"/></p>

</div>
<div id='tab-get-signatures-from-hierarchical-partition-3'>
<pre><code class="r">get_signatures(res_rh, merge_node = merge_node_param(min_n_signatures = 2054))
</code></pre>

<p><img src="figure_cola/tab-get-signatures-from-hierarchical-partition-3-1.png" alt="plot of chunk tab-get-signatures-from-hierarchical-partition-3"/></p>

</div>
</div>




Compare signatures from different nodes:


```r
compare_signatures(res_rh, verbose = FALSE)
```

![plot of chunk unnamed-chunk-24](figure_cola/unnamed-chunk-24-1.png)

If there are too many signatures, `top_signatures = ...` can be set to only show the 
signatures with the highest FDRs. Note it only works on every node and the final signatures
are the union of all signatures of all nodes.


```r
# code only for demonstration
# e.g. to show the top 500 most significant rows on each node.
tb = get_signature(res_rh, top_signatures = 500)
```


### Test to known annotations

Test correlation between subgroups and known annotations. If the known
annotation is numeric, one-way ANOVA test is applied, and if the known
annotation is discrete, chi-squared contingency table test is applied.




<script>
$( function() {
	$( '#tabs-test-to-known-factors-from-hierarchical-partition' ).tabs();
} );
</script>
<div id='tabs-test-to-known-factors-from-hierarchical-partition'>
<ul>
<li><a href='#tab-test-to-known-factors-from-hierarchical-partition-1'>n_signatures ≥ 134</a></li>
<li><a href='#tab-test-to-known-factors-from-hierarchical-partition-2'>n_signatures ≥ 650</a></li>
<li><a href='#tab-test-to-known-factors-from-hierarchical-partition-3'>n_signatures ≥ 2054</a></li>
</ul>
<div id='tab-test-to-known-factors-from-hierarchical-partition-1'>
<pre><code class="r">test_to_known_factors(res_rh, merge_node = merge_node_param(min_n_signatures = 134))
</code></pre>

<pre><code>#&gt;        ALL.AML
#&gt; class 5.34e-09
</code></pre>

</div>
<div id='tab-test-to-known-factors-from-hierarchical-partition-2'>
<pre><code class="r">test_to_known_factors(res_rh, merge_node = merge_node_param(min_n_signatures = 650))
</code></pre>

<pre><code>#&gt;       ALL.AML
#&gt; class 2.8e-09
</code></pre>

</div>
<div id='tab-test-to-known-factors-from-hierarchical-partition-3'>
<pre><code class="r">test_to_known_factors(res_rh, merge_node = merge_node_param(min_n_signatures = 2054))
</code></pre>

<pre><code>#&gt;        ALL.AML
#&gt; class 1.88e-10
</code></pre>

</div>
</div>



## Results for each node


---------------------------------------------------




### Node0


Child nodes: 
                [Node01](#Node01)
        ,
                [Node02](#Node02)
        ,
                Node03-leaf
        .







The object with results only for a single top-value method and a single partitioning method 
can be extracted as:

```r
res = res_rh["0"]
```

A summary of `res` and all the functions that can be applied to it:

```r
res
```

```
#> A 'ConsensusPartition' object with k = 2, 3, 4.
#>   On a matrix with 3910 rows and 72 columns.
#>   Top rows (391) are extracted by 'ATC' method.
#>   Subgroups are detected by 'skmeans' method.
#>   Performed in total 150 partitions by row resampling.
#>   Best k for subgroups seems to be 3.
#> 
#> Following methods can be applied to this 'ConsensusPartition' object:
#>  [1] "cola_report"             "collect_classes"         "collect_plots"          
#>  [4] "collect_stats"           "colnames"                "compare_partitions"     
#>  [7] "compare_signatures"      "consensus_heatmap"       "dimension_reduction"    
#> [10] "functional_enrichment"   "get_anno_col"            "get_anno"               
#> [13] "get_classes"             "get_consensus"           "get_matrix"             
#> [16] "get_membership"          "get_param"               "get_signatures"         
#> [19] "get_stats"               "is_best_k"               "is_stable_k"            
#> [22] "membership_heatmap"      "ncol"                    "nrow"                   
#> [25] "plot_ecdf"               "predict_classes"         "rownames"               
#> [28] "select_partition_number" "show"                    "suggest_best_k"         
#> [31] "test_to_known_factors"   "top_rows_heatmap"
```

`collect_plots()` function collects all the plots made from `res` for all `k` (number of subgroups)
into one single page to provide an easy and fast comparison between different `k`.

```r
collect_plots(res)
```

![plot of chunk node-0-collect-plots](figure_cola/node-0-collect-plots-1.png)

The plots are:

- The first row: a plot of the eCDF (empirical cumulative distribution
  function) curves of the consensus matrix for each `k` and the heatmap of
  predicted classes for each `k`.
- The second row: heatmaps of the consensus matrix for each `k`.
- The third row: heatmaps of the membership matrix for each `k`.
- The fouth row: heatmaps of the signatures for each `k`.

All the plots in panels can be made by individual functions and they are
plotted later in this section.

`select_partition_number()` produces several plots showing different
statistics for choosing "optimized" `k`. There are following statistics:

- eCDF curves of the consensus matrix for each `k`;
- 1-PAC. [The PAC score](https://en.wikipedia.org/wiki/Consensus_clustering#Over-interpretation_potential_of_consensus_clustering)
  measures the proportion of the ambiguous subgrouping.
- Mean silhouette score.
- Concordance. The mean probability of fiting the consensus subgroup labels in all
  partitions.
- Area increased. Denote $A_k$ as the area under the eCDF curve for current
  `k`, the area increased is defined as $A_k - A_{k-1}$.
- Rand index. The percent of pairs of samples that are both in a same cluster
  or both are not in a same cluster in the partition of k and k-1.
- Jaccard index. The ratio of pairs of samples are both in a same cluster in
  the partition of k and k-1 and the pairs of samples are both in a same
  cluster in the partition k or k-1.

The detailed explanations of these statistics can be found in [the _cola_
vignette](https://jokergoo.github.io/cola_vignettes/cola.html#toc_13).

Generally speaking, higher 1-PAC score, higher mean silhouette score or higher
concordance corresponds to better partition. Rand index and Jaccard index
measure how similar the current partition is compared to partition with `k-1`.
If they are too similar, we won't accept `k` is better than `k-1`.

```r
select_partition_number(res)
```

![plot of chunk node-0-select-partition-number](figure_cola/node-0-select-partition-number-1.png)

The numeric values for all these statistics can be obtained by `get_stats()`.

```r
get_stats(res)
```

```
#>   k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#> 2 2 1.000           0.980       0.992          0.507 0.493   0.493
#> 3 3 1.000           0.968       0.984          0.251 0.869   0.737
#> 4 4 0.877           0.888       0.915          0.111 0.915   0.772
```

`suggest_best_k()` suggests the best $k$ based on these statistics. The rules are as follows:

- All $k$ with Jaccard index larger than 0.95 are removed because increasing
  $k$ does not provide enough extra information. If all $k$ are removed, it is
  marked as no subgroup is detected.
- For all $k$ with 1-PAC score larger than 0.9, the maximal $k$ is taken as
  the best $k$, and other $k$ are marked as optional $k$.
- If it does not fit the second rule. The $k$ with the maximal vote of the
  highest 1-PAC score, highest mean silhouette, and highest concordance is
  taken as the best $k$.

```r
suggest_best_k(res)
```

```
#> [1] 3
#> attr(,"optional")
#> [1] 2
```

There is also optional best $k$ = 2 that is worth to check.

Following is the table of the partitions (You need to click the **show/hide
code output** link to see it). The membership matrix (columns with name `p*`)
is inferred by
[`clue::cl_consensus()`](https://www.rdocumentation.org/link/cl_consensus?package=clue)
function with the `SE` method. Basically the value in the membership matrix
represents the probability to belong to a certain group. The finall subgroup
label for an item is determined with the group with highest probability it
belongs to.

In `get_classes()` function, the entropy is calculated from the membership
matrix and the silhouette score is calculated from the consensus matrix.



<script>
$( function() {
	$( '#tabs-node-0-get-classes' ).tabs();
} );
</script>
<div id='tabs-node-0-get-classes'>
<ul>
<li><a href='#tab-node-0-get-classes-1'>k = 2</a></li>
<li><a href='#tab-node-0-get-classes-2'>k = 3</a></li>
<li><a href='#tab-node-0-get-classes-3'>k = 4</a></li>
</ul>

<div id='tab-node-0-get-classes-1'>
<p><a id='tab-node-0-get-classes-1-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 2), get_membership(res, k = 2))
</code></pre>

<pre><code>#&gt;           class entropy silhouette   p1   p2
#&gt; sample_39     2   0.000      0.999 0.00 1.00
#&gt; sample_40     2   0.000      0.999 0.00 1.00
#&gt; sample_42     1   0.000      0.984 1.00 0.00
#&gt; sample_47     1   0.000      0.984 1.00 0.00
#&gt; sample_48     1   0.000      0.984 1.00 0.00
#&gt; sample_49     2   0.000      0.999 0.00 1.00
#&gt; sample_41     1   0.000      0.984 1.00 0.00
#&gt; sample_43     1   0.000      0.984 1.00 0.00
#&gt; sample_44     1   0.000      0.984 1.00 0.00
#&gt; sample_45     1   0.000      0.984 1.00 0.00
#&gt; sample_46     1   0.000      0.984 1.00 0.00
#&gt; sample_70     1   0.000      0.984 1.00 0.00
#&gt; sample_71     1   0.000      0.984 1.00 0.00
#&gt; sample_72     1   0.000      0.984 1.00 0.00
#&gt; sample_68     1   0.000      0.984 1.00 0.00
#&gt; sample_69     1   0.000      0.984 1.00 0.00
#&gt; sample_67     1   0.000      0.984 1.00 0.00
#&gt; sample_55     2   0.000      0.999 0.00 1.00
#&gt; sample_56     2   0.000      0.999 0.00 1.00
#&gt; sample_59     1   0.000      0.984 1.00 0.00
#&gt; sample_52     2   0.000      0.999 0.00 1.00
#&gt; sample_53     2   0.000      0.999 0.00 1.00
#&gt; sample_51     2   0.000      0.999 0.00 1.00
#&gt; sample_50     2   0.000      0.999 0.00 1.00
#&gt; sample_54     1   0.000      0.984 1.00 0.00
#&gt; sample_57     2   0.000      0.999 0.00 1.00
#&gt; sample_58     2   0.000      0.999 0.00 1.00
#&gt; sample_60     1   0.584      0.833 0.86 0.14
#&gt; sample_61     2   0.000      0.999 0.00 1.00
#&gt; sample_65     2   0.000      0.999 0.00 1.00
#&gt; sample_66     1   0.000      0.984 1.00 0.00
#&gt; sample_63     2   0.000      0.999 0.00 1.00
#&gt; sample_64     2   0.000      0.999 0.00 1.00
#&gt; sample_62     2   0.000      0.999 0.00 1.00
#&gt; sample_1      2   0.242      0.957 0.04 0.96
#&gt; sample_2      1   0.000      0.984 1.00 0.00
#&gt; sample_3      2   0.000      0.999 0.00 1.00
#&gt; sample_4      1   0.971      0.339 0.60 0.40
#&gt; sample_5      1   0.000      0.984 1.00 0.00
#&gt; sample_6      2   0.000      0.999 0.00 1.00
#&gt; sample_7      2   0.000      0.999 0.00 1.00
#&gt; sample_8      2   0.000      0.999 0.00 1.00
#&gt; sample_9      1   0.000      0.984 1.00 0.00
#&gt; sample_10     1   0.000      0.984 1.00 0.00
#&gt; sample_11     1   0.000      0.984 1.00 0.00
#&gt; sample_12     2   0.000      0.999 0.00 1.00
#&gt; sample_13     1   0.000      0.984 1.00 0.00
#&gt; sample_14     1   0.000      0.984 1.00 0.00
#&gt; sample_15     1   0.000      0.984 1.00 0.00
#&gt; sample_16     1   0.000      0.984 1.00 0.00
#&gt; sample_17     1   0.000      0.984 1.00 0.00
#&gt; sample_18     1   0.000      0.984 1.00 0.00
#&gt; sample_19     1   0.000      0.984 1.00 0.00
#&gt; sample_20     1   0.000      0.984 1.00 0.00
#&gt; sample_21     1   0.000      0.984 1.00 0.00
#&gt; sample_22     2   0.000      0.999 0.00 1.00
#&gt; sample_23     2   0.000      0.999 0.00 1.00
#&gt; sample_24     1   0.000      0.984 1.00 0.00
#&gt; sample_25     2   0.000      0.999 0.00 1.00
#&gt; sample_26     1   0.000      0.984 1.00 0.00
#&gt; sample_27     2   0.000      0.999 0.00 1.00
#&gt; sample_34     2   0.000      0.999 0.00 1.00
#&gt; sample_35     2   0.000      0.999 0.00 1.00
#&gt; sample_36     2   0.000      0.999 0.00 1.00
#&gt; sample_37     2   0.000      0.999 0.00 1.00
#&gt; sample_38     2   0.000      0.999 0.00 1.00
#&gt; sample_28     2   0.000      0.999 0.00 1.00
#&gt; sample_29     1   0.000      0.984 1.00 0.00
#&gt; sample_30     2   0.000      0.999 0.00 1.00
#&gt; sample_31     2   0.000      0.999 0.00 1.00
#&gt; sample_32     2   0.000      0.999 0.00 1.00
#&gt; sample_33     2   0.000      0.999 0.00 1.00
</code></pre>

<script>
$('#tab-node-0-get-classes-1-a').parent().next().next().hide();
$('#tab-node-0-get-classes-1-a').click(function(){
  $('#tab-node-0-get-classes-1-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-0-get-classes-2'>
<p><a id='tab-node-0-get-classes-2-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 3), get_membership(res, k = 3))
</code></pre>

<pre><code>#&gt;           class entropy silhouette   p1   p2   p3
#&gt; sample_39     3  0.2066      0.951 0.00 0.06 0.94
#&gt; sample_40     3  0.0892      0.989 0.00 0.02 0.98
#&gt; sample_42     1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_47     1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_48     1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_49     3  0.0892      0.989 0.00 0.02 0.98
#&gt; sample_41     1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_43     1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_44     1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_45     1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_46     1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_70     1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_71     1  0.2537      0.894 0.92 0.08 0.00
#&gt; sample_72     1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_68     1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_69     1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_67     1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_55     3  0.0000      0.981 0.00 0.00 1.00
#&gt; sample_56     3  0.0892      0.989 0.00 0.02 0.98
#&gt; sample_59     1  0.2537      0.903 0.92 0.00 0.08
#&gt; sample_52     2  0.0000      1.000 0.00 1.00 0.00
#&gt; sample_53     2  0.0000      1.000 0.00 1.00 0.00
#&gt; sample_51     2  0.0000      1.000 0.00 1.00 0.00
#&gt; sample_50     2  0.0000      1.000 0.00 1.00 0.00
#&gt; sample_54     1  0.0892      0.962 0.98 0.00 0.02
#&gt; sample_57     2  0.0000      1.000 0.00 1.00 0.00
#&gt; sample_58     2  0.0000      1.000 0.00 1.00 0.00
#&gt; sample_60     1  0.8143      0.322 0.56 0.08 0.36
#&gt; sample_61     2  0.0000      1.000 0.00 1.00 0.00
#&gt; sample_65     2  0.0000      1.000 0.00 1.00 0.00
#&gt; sample_66     1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_63     2  0.0000      1.000 0.00 1.00 0.00
#&gt; sample_64     2  0.0000      1.000 0.00 1.00 0.00
#&gt; sample_62     2  0.0000      1.000 0.00 1.00 0.00
#&gt; sample_1      3  0.0000      0.981 0.00 0.00 1.00
#&gt; sample_2      1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_3      3  0.0000      0.981 0.00 0.00 1.00
#&gt; sample_4      3  0.0000      0.981 0.00 0.00 1.00
#&gt; sample_5      1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_6      3  0.0892      0.989 0.00 0.02 0.98
#&gt; sample_7      3  0.0892      0.989 0.00 0.02 0.98
#&gt; sample_8      3  0.0892      0.989 0.00 0.02 0.98
#&gt; sample_9      1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_10     1  0.0892      0.962 0.98 0.00 0.02
#&gt; sample_11     1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_12     2  0.0000      1.000 0.00 1.00 0.00
#&gt; sample_13     1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_14     1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_15     1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_16     1  0.0892      0.962 0.98 0.00 0.02
#&gt; sample_17     1  0.0892      0.962 0.98 0.00 0.02
#&gt; sample_18     1  0.4796      0.738 0.78 0.00 0.22
#&gt; sample_19     1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_20     1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_21     1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_22     2  0.0000      1.000 0.00 1.00 0.00
#&gt; sample_23     3  0.0892      0.989 0.00 0.02 0.98
#&gt; sample_24     1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_25     2  0.0000      1.000 0.00 1.00 0.00
#&gt; sample_26     1  0.2066      0.930 0.94 0.00 0.06
#&gt; sample_27     3  0.0892      0.989 0.00 0.02 0.98
#&gt; sample_34     2  0.0000      1.000 0.00 1.00 0.00
#&gt; sample_35     2  0.0000      1.000 0.00 1.00 0.00
#&gt; sample_36     2  0.0000      1.000 0.00 1.00 0.00
#&gt; sample_37     2  0.0000      1.000 0.00 1.00 0.00
#&gt; sample_38     2  0.0000      1.000 0.00 1.00 0.00
#&gt; sample_28     2  0.0000      1.000 0.00 1.00 0.00
#&gt; sample_29     1  0.0000      0.972 1.00 0.00 0.00
#&gt; sample_30     2  0.0000      1.000 0.00 1.00 0.00
#&gt; sample_31     2  0.0000      1.000 0.00 1.00 0.00
#&gt; sample_32     2  0.0000      1.000 0.00 1.00 0.00
#&gt; sample_33     2  0.0000      1.000 0.00 1.00 0.00
</code></pre>

<script>
$('#tab-node-0-get-classes-2-a').parent().next().next().hide();
$('#tab-node-0-get-classes-2-a').click(function(){
  $('#tab-node-0-get-classes-2-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-0-get-classes-3'>
<p><a id='tab-node-0-get-classes-3-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 4), get_membership(res, k = 4))
</code></pre>

<pre><code>#&gt;           class entropy silhouette   p1   p2   p3   p4
#&gt; sample_39     3  0.0000      0.987 0.00 0.00 1.00 0.00
#&gt; sample_40     3  0.0000      0.987 0.00 0.00 1.00 0.00
#&gt; sample_42     4  0.4948      0.932 0.44 0.00 0.00 0.56
#&gt; sample_47     1  0.0000      0.846 1.00 0.00 0.00 0.00
#&gt; sample_48     1  0.0000      0.846 1.00 0.00 0.00 0.00
#&gt; sample_49     3  0.0000      0.987 0.00 0.00 1.00 0.00
#&gt; sample_41     1  0.0000      0.846 1.00 0.00 0.00 0.00
#&gt; sample_43     1  0.1211      0.812 0.96 0.00 0.00 0.04
#&gt; sample_44     1  0.0000      0.846 1.00 0.00 0.00 0.00
#&gt; sample_45     1  0.0000      0.846 1.00 0.00 0.00 0.00
#&gt; sample_46     1  0.1211      0.807 0.96 0.00 0.00 0.04
#&gt; sample_70     1  0.3400      0.661 0.82 0.00 0.00 0.18
#&gt; sample_71     4  0.5428      0.891 0.38 0.02 0.00 0.60
#&gt; sample_72     4  0.4907      0.931 0.42 0.00 0.00 0.58
#&gt; sample_68     1  0.0000      0.846 1.00 0.00 0.00 0.00
#&gt; sample_69     1  0.0000      0.846 1.00 0.00 0.00 0.00
#&gt; sample_67     4  0.4907      0.931 0.42 0.00 0.00 0.58
#&gt; sample_55     3  0.2647      0.890 0.00 0.00 0.88 0.12
#&gt; sample_56     3  0.0000      0.987 0.00 0.00 1.00 0.00
#&gt; sample_59     1  0.4948      0.371 0.56 0.00 0.00 0.44
#&gt; sample_52     2  0.0000      0.996 0.00 1.00 0.00 0.00
#&gt; sample_53     2  0.0000      0.996 0.00 1.00 0.00 0.00
#&gt; sample_51     2  0.0000      0.996 0.00 1.00 0.00 0.00
#&gt; sample_50     2  0.0000      0.996 0.00 1.00 0.00 0.00
#&gt; sample_54     1  0.4790      0.435 0.62 0.00 0.00 0.38
#&gt; sample_57     2  0.0000      0.996 0.00 1.00 0.00 0.00
#&gt; sample_58     2  0.0000      0.996 0.00 1.00 0.00 0.00
#&gt; sample_60     1  0.7192      0.248 0.48 0.02 0.08 0.42
#&gt; sample_61     2  0.0000      0.996 0.00 1.00 0.00 0.00
#&gt; sample_65     2  0.0000      0.996 0.00 1.00 0.00 0.00
#&gt; sample_66     1  0.0000      0.846 1.00 0.00 0.00 0.00
#&gt; sample_63     2  0.0000      0.996 0.00 1.00 0.00 0.00
#&gt; sample_64     2  0.0000      0.996 0.00 1.00 0.00 0.00
#&gt; sample_62     2  0.0000      0.996 0.00 1.00 0.00 0.00
#&gt; sample_1      3  0.0000      0.987 0.00 0.00 1.00 0.00
#&gt; sample_2      4  0.4977      0.918 0.46 0.00 0.00 0.54
#&gt; sample_3      3  0.0707      0.979 0.00 0.00 0.98 0.02
#&gt; sample_4      3  0.0000      0.987 0.00 0.00 1.00 0.00
#&gt; sample_5      1  0.0000      0.846 1.00 0.00 0.00 0.00
#&gt; sample_6      3  0.0707      0.979 0.00 0.00 0.98 0.02
#&gt; sample_7      3  0.0000      0.987 0.00 0.00 1.00 0.00
#&gt; sample_8      3  0.0000      0.987 0.00 0.00 1.00 0.00
#&gt; sample_9      1  0.0000      0.846 1.00 0.00 0.00 0.00
#&gt; sample_10     4  0.4713      0.816 0.36 0.00 0.00 0.64
#&gt; sample_11     4  0.4994      0.889 0.48 0.00 0.00 0.52
#&gt; sample_12     2  0.0000      0.996 0.00 1.00 0.00 0.00
#&gt; sample_13     1  0.0000      0.846 1.00 0.00 0.00 0.00
#&gt; sample_14     1  0.0000      0.846 1.00 0.00 0.00 0.00
#&gt; sample_15     1  0.0000      0.846 1.00 0.00 0.00 0.00
#&gt; sample_16     1  0.0000      0.846 1.00 0.00 0.00 0.00
#&gt; sample_17     1  0.0707      0.828 0.98 0.00 0.00 0.02
#&gt; sample_18     1  0.4948      0.352 0.56 0.00 0.00 0.44
#&gt; sample_19     1  0.0000      0.846 1.00 0.00 0.00 0.00
#&gt; sample_20     1  0.0000      0.846 1.00 0.00 0.00 0.00
#&gt; sample_21     1  0.0000      0.846 1.00 0.00 0.00 0.00
#&gt; sample_22     2  0.0000      0.996 0.00 1.00 0.00 0.00
#&gt; sample_23     3  0.0707      0.979 0.00 0.00 0.98 0.02
#&gt; sample_24     1  0.0000      0.846 1.00 0.00 0.00 0.00
#&gt; sample_25     2  0.2345      0.889 0.00 0.90 0.10 0.00
#&gt; sample_26     1  0.5271      0.456 0.64 0.00 0.02 0.34
#&gt; sample_27     3  0.0000      0.987 0.00 0.00 1.00 0.00
#&gt; sample_34     2  0.0000      0.996 0.00 1.00 0.00 0.00
#&gt; sample_35     2  0.0000      0.996 0.00 1.00 0.00 0.00
#&gt; sample_36     2  0.0000      0.996 0.00 1.00 0.00 0.00
#&gt; sample_37     2  0.0000      0.996 0.00 1.00 0.00 0.00
#&gt; sample_38     2  0.0000      0.996 0.00 1.00 0.00 0.00
#&gt; sample_28     2  0.0000      0.996 0.00 1.00 0.00 0.00
#&gt; sample_29     4  0.4948      0.932 0.44 0.00 0.00 0.56
#&gt; sample_30     2  0.0000      0.996 0.00 1.00 0.00 0.00
#&gt; sample_31     2  0.0000      0.996 0.00 1.00 0.00 0.00
#&gt; sample_32     2  0.0000      0.996 0.00 1.00 0.00 0.00
#&gt; sample_33     2  0.0000      0.996 0.00 1.00 0.00 0.00
</code></pre>

<script>
$('#tab-node-0-get-classes-3-a').parent().next().next().hide();
$('#tab-node-0-get-classes-3-a').click(function(){
  $('#tab-node-0-get-classes-3-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>
</div>

Heatmaps for the consensus matrix. It visualizes the probability of two
samples to be in a same group.




<script>
$( function() {
	$( '#tabs-node-0-consensus-heatmap' ).tabs();
} );
</script>
<div id='tabs-node-0-consensus-heatmap'>
<ul>
<li><a href='#tab-node-0-consensus-heatmap-1'>k = 2</a></li>
<li><a href='#tab-node-0-consensus-heatmap-2'>k = 3</a></li>
<li><a href='#tab-node-0-consensus-heatmap-3'>k = 4</a></li>
</ul>
<div id='tab-node-0-consensus-heatmap-1'>
<pre><code class="r">consensus_heatmap(res, k = 2)
</code></pre>

<p><img src="figure_cola/tab-node-0-consensus-heatmap-1-1.png" alt="plot of chunk tab-node-0-consensus-heatmap-1"/></p>

</div>
<div id='tab-node-0-consensus-heatmap-2'>
<pre><code class="r">consensus_heatmap(res, k = 3)
</code></pre>

<p><img src="figure_cola/tab-node-0-consensus-heatmap-2-1.png" alt="plot of chunk tab-node-0-consensus-heatmap-2"/></p>

</div>
<div id='tab-node-0-consensus-heatmap-3'>
<pre><code class="r">consensus_heatmap(res, k = 4)
</code></pre>

<p><img src="figure_cola/tab-node-0-consensus-heatmap-3-1.png" alt="plot of chunk tab-node-0-consensus-heatmap-3"/></p>

</div>
</div>

Heatmaps for the membership of samples in all partitions to see how consistent they are:





<script>
$( function() {
	$( '#tabs-node-0-membership-heatmap' ).tabs();
} );
</script>
<div id='tabs-node-0-membership-heatmap'>
<ul>
<li><a href='#tab-node-0-membership-heatmap-1'>k = 2</a></li>
<li><a href='#tab-node-0-membership-heatmap-2'>k = 3</a></li>
<li><a href='#tab-node-0-membership-heatmap-3'>k = 4</a></li>
</ul>
<div id='tab-node-0-membership-heatmap-1'>
<pre><code class="r">membership_heatmap(res, k = 2)
</code></pre>

<p><img src="figure_cola/tab-node-0-membership-heatmap-1-1.png" alt="plot of chunk tab-node-0-membership-heatmap-1"/></p>

</div>
<div id='tab-node-0-membership-heatmap-2'>
<pre><code class="r">membership_heatmap(res, k = 3)
</code></pre>

<p><img src="figure_cola/tab-node-0-membership-heatmap-2-1.png" alt="plot of chunk tab-node-0-membership-heatmap-2"/></p>

</div>
<div id='tab-node-0-membership-heatmap-3'>
<pre><code class="r">membership_heatmap(res, k = 4)
</code></pre>

<p><img src="figure_cola/tab-node-0-membership-heatmap-3-1.png" alt="plot of chunk tab-node-0-membership-heatmap-3"/></p>

</div>
</div>

As soon as the classes for columns are determined, the signatures
that are significantly different between subgroups can be looked for. 
Following are the heatmaps for signatures.




Signature heatmaps where rows are scaled:



<script>
$( function() {
	$( '#tabs-node-0-get-signatures' ).tabs();
} );
</script>
<div id='tabs-node-0-get-signatures'>
<ul>
<li><a href='#tab-node-0-get-signatures-1'>k = 2</a></li>
<li><a href='#tab-node-0-get-signatures-2'>k = 3</a></li>
<li><a href='#tab-node-0-get-signatures-3'>k = 4</a></li>
</ul>
<div id='tab-node-0-get-signatures-1'>
<pre><code class="r">get_signatures(res, k = 2)
</code></pre>

<p><img src="figure_cola/tab-node-0-get-signatures-1-1.png" alt="plot of chunk tab-node-0-get-signatures-1"/></p>

</div>
<div id='tab-node-0-get-signatures-2'>
<pre><code class="r">get_signatures(res, k = 3)
</code></pre>

<p><img src="figure_cola/tab-node-0-get-signatures-2-1.png" alt="plot of chunk tab-node-0-get-signatures-2"/></p>

</div>
<div id='tab-node-0-get-signatures-3'>
<pre><code class="r">get_signatures(res, k = 4)
</code></pre>

<p><img src="figure_cola/tab-node-0-get-signatures-3-1.png" alt="plot of chunk tab-node-0-get-signatures-3"/></p>

</div>
</div>



Signature heatmaps where rows are not scaled:


<script>
$( function() {
	$( '#tabs-node-0-get-signatures-no-scale' ).tabs();
} );
</script>
<div id='tabs-node-0-get-signatures-no-scale'>
<ul>
<li><a href='#tab-node-0-get-signatures-no-scale-1'>k = 2</a></li>
<li><a href='#tab-node-0-get-signatures-no-scale-2'>k = 3</a></li>
<li><a href='#tab-node-0-get-signatures-no-scale-3'>k = 4</a></li>
</ul>
<div id='tab-node-0-get-signatures-no-scale-1'>
<pre><code class="r">get_signatures(res, k = 2, scale_rows = FALSE)
</code></pre>

<p><img src="figure_cola/tab-node-0-get-signatures-no-scale-1-1.png" alt="plot of chunk tab-node-0-get-signatures-no-scale-1"/></p>

</div>
<div id='tab-node-0-get-signatures-no-scale-2'>
<pre><code class="r">get_signatures(res, k = 3, scale_rows = FALSE)
</code></pre>

<p><img src="figure_cola/tab-node-0-get-signatures-no-scale-2-1.png" alt="plot of chunk tab-node-0-get-signatures-no-scale-2"/></p>

</div>
<div id='tab-node-0-get-signatures-no-scale-3'>
<pre><code class="r">get_signatures(res, k = 4, scale_rows = FALSE)
</code></pre>

<p><img src="figure_cola/tab-node-0-get-signatures-no-scale-3-1.png" alt="plot of chunk tab-node-0-get-signatures-no-scale-3"/></p>

</div>
</div>



Compare the overlap of signatures from different k:

```r
compare_signatures(res)
```

![plot of chunk node-0-signature_compare](figure_cola/node-0-signature_compare-1.png)

`get_signature()` returns a data frame invisibly. To get the list of signatures, the function
call should be assigned to a variable explicitly. In following code, if `plot` argument is set
to `FALSE`, no heatmap is plotted while only the differential analysis is performed.

```r
# code only for demonstration
tb = get_signature(res, k = ..., plot = FALSE)
```

An example of the output of `tb` is:

```
#>   which_row         fdr    mean_1    mean_2 scaled_mean_1 scaled_mean_2 km
#> 1        38 0.042760348  8.373488  9.131774    -0.5533452     0.5164555  1
#> 2        40 0.018707592  7.106213  8.469186    -0.6173731     0.5762149  1
#> 3        55 0.019134737 10.221463 11.207825    -0.6159697     0.5749050  1
#> 4        59 0.006059896  5.921854  7.869574    -0.6899429     0.6439467  1
#> 5        60 0.018055526  8.928898 10.211722    -0.6204761     0.5791110  1
#> 6        98 0.009384629 15.714769 14.887706     0.6635654    -0.6193277  2
...
```

The columns in `tb` are:

1. `which_row`: row indices corresponding to the input matrix.
2. `fdr`: FDR for the differential test. 
3. `mean_x`: The mean value in group x.
4. `scaled_mean_x`: The mean value in group x after rows are scaled.
5. `km`: Row groups if k-means clustering is applied to rows (which is done by automatically selecting number of clusters).

If there are too many signatures, `top_signatures = ...` can be set to only show the 
signatures with the highest FDRs:

```r
# code only for demonstration
# e.g. to show the top 500 most significant rows
tb = get_signature(res, k = ..., top_signatures = 500)
```

If the signatures are defined as these which are uniquely high in current group, `diff_method` argument
can be set to `"uniquely_high_in_one_group"`:

```r
# code only for demonstration
tb = get_signature(res, k = ..., diff_method = "uniquely_high_in_one_group")
```




UMAP plot which shows how samples are separated.


<script>
$( function() {
	$( '#tabs-node-0-dimension-reduction' ).tabs();
} );
</script>
<div id='tabs-node-0-dimension-reduction'>
<ul>
<li><a href='#tab-node-0-dimension-reduction-1'>k = 2</a></li>
<li><a href='#tab-node-0-dimension-reduction-2'>k = 3</a></li>
<li><a href='#tab-node-0-dimension-reduction-3'>k = 4</a></li>
</ul>
<div id='tab-node-0-dimension-reduction-1'>
<pre><code class="r">dimension_reduction(res, k = 2, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-0-dimension-reduction-1-1.png" alt="plot of chunk tab-node-0-dimension-reduction-1"/></p>

</div>
<div id='tab-node-0-dimension-reduction-2'>
<pre><code class="r">dimension_reduction(res, k = 3, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-0-dimension-reduction-2-1.png" alt="plot of chunk tab-node-0-dimension-reduction-2"/></p>

</div>
<div id='tab-node-0-dimension-reduction-3'>
<pre><code class="r">dimension_reduction(res, k = 4, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-0-dimension-reduction-3-1.png" alt="plot of chunk tab-node-0-dimension-reduction-3"/></p>

</div>
</div>



Following heatmap shows how subgroups are split when increasing `k`:

```r
collect_classes(res)
```

![plot of chunk node-0-collect-classes](figure_cola/node-0-collect-classes-1.png)




Test correlation between subgroups and known annotations. If the known
annotation is numeric, one-way ANOVA test is applied, and if the known
annotation is discrete, chi-squared contingency table test is applied.

```r
test_to_known_factors(res)
```

```
#>             n_sample ALL.AML(p-value) k
#> ATC:skmeans       71         1.01e-04 2
#> ATC:skmeans       71         6.09e-11 3
#> ATC:skmeans       67         3.09e-10 4
```




If matrix rows can be associated to genes, consider to use `functional_enrichment(res,
...)` to perform function enrichment for the signature genes. See [this vignette](https://jokergoo.github.io/cola_vignettes/functional_enrichment.html) for more detailed explanations.


 

---------------------------------------------------




### Node01


Parent node: [Node0](#Node0).
Child nodes: 
                Node011-leaf
        ,
                Node012-leaf
        ,
                Node013-leaf
        ,
                Node021-leaf
        ,
                Node022-leaf
        .







The object with results only for a single top-value method and a single partitioning method 
can be extracted as:

```r
res = res_rh["01"]
```

A summary of `res` and all the functions that can be applied to it:

```r
res
```

```
#> A 'ConsensusPartition' object with k = 2, 3, 4.
#>   On a matrix with 3910 rows and 35 columns.
#>   Top rows (391) are extracted by 'ATC' method.
#>   Subgroups are detected by 'skmeans' method.
#>   Performed in total 150 partitions by row resampling.
#>   Best k for subgroups seems to be 3.
#> 
#> Following methods can be applied to this 'ConsensusPartition' object:
#>  [1] "cola_report"             "collect_classes"         "collect_plots"          
#>  [4] "collect_stats"           "colnames"                "compare_partitions"     
#>  [7] "compare_signatures"      "consensus_heatmap"       "dimension_reduction"    
#> [10] "functional_enrichment"   "get_anno_col"            "get_anno"               
#> [13] "get_classes"             "get_consensus"           "get_matrix"             
#> [16] "get_membership"          "get_param"               "get_signatures"         
#> [19] "get_stats"               "is_best_k"               "is_stable_k"            
#> [22] "membership_heatmap"      "ncol"                    "nrow"                   
#> [25] "plot_ecdf"               "predict_classes"         "rownames"               
#> [28] "select_partition_number" "show"                    "suggest_best_k"         
#> [31] "test_to_known_factors"   "top_rows_heatmap"
```

`collect_plots()` function collects all the plots made from `res` for all `k` (number of subgroups)
into one single page to provide an easy and fast comparison between different `k`.

```r
collect_plots(res)
```

![plot of chunk node-01-collect-plots](figure_cola/node-01-collect-plots-1.png)

The plots are:

- The first row: a plot of the eCDF (empirical cumulative distribution
  function) curves of the consensus matrix for each `k` and the heatmap of
  predicted classes for each `k`.
- The second row: heatmaps of the consensus matrix for each `k`.
- The third row: heatmaps of the membership matrix for each `k`.
- The fouth row: heatmaps of the signatures for each `k`.

All the plots in panels can be made by individual functions and they are
plotted later in this section.

`select_partition_number()` produces several plots showing different
statistics for choosing "optimized" `k`. There are following statistics:

- eCDF curves of the consensus matrix for each `k`;
- 1-PAC. [The PAC score](https://en.wikipedia.org/wiki/Consensus_clustering#Over-interpretation_potential_of_consensus_clustering)
  measures the proportion of the ambiguous subgrouping.
- Mean silhouette score.
- Concordance. The mean probability of fiting the consensus subgroup labels in all
  partitions.
- Area increased. Denote $A_k$ as the area under the eCDF curve for current
  `k`, the area increased is defined as $A_k - A_{k-1}$.
- Rand index. The percent of pairs of samples that are both in a same cluster
  or both are not in a same cluster in the partition of k and k-1.
- Jaccard index. The ratio of pairs of samples are both in a same cluster in
  the partition of k and k-1 and the pairs of samples are both in a same
  cluster in the partition k or k-1.

The detailed explanations of these statistics can be found in [the _cola_
vignette](https://jokergoo.github.io/cola_vignettes/cola.html#toc_13).

Generally speaking, higher 1-PAC score, higher mean silhouette score or higher
concordance corresponds to better partition. Rand index and Jaccard index
measure how similar the current partition is compared to partition with `k-1`.
If they are too similar, we won't accept `k` is better than `k-1`.

```r
select_partition_number(res)
```

![plot of chunk node-01-select-partition-number](figure_cola/node-01-select-partition-number-1.png)

The numeric values for all these statistics can be obtained by `get_stats()`.

```r
get_stats(res)
```

```
#>   k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#> 2 2 1.000           0.973       0.988         0.4982 0.506   0.506
#> 3 3 1.000           0.954       0.983         0.3745 0.703   0.468
#> 4 4 0.805           0.767       0.897         0.0993 0.904   0.708
```

`suggest_best_k()` suggests the best $k$ based on these statistics. The rules are as follows:

- All $k$ with Jaccard index larger than 0.95 are removed because increasing
  $k$ does not provide enough extra information. If all $k$ are removed, it is
  marked as no subgroup is detected.
- For all $k$ with 1-PAC score larger than 0.9, the maximal $k$ is taken as
  the best $k$, and other $k$ are marked as optional $k$.
- If it does not fit the second rule. The $k$ with the maximal vote of the
  highest 1-PAC score, highest mean silhouette, and highest concordance is
  taken as the best $k$.

```r
suggest_best_k(res)
```

```
#> [1] 3
#> attr(,"optional")
#> [1] 2
```

There is also optional best $k$ = 2 that is worth to check.

Following is the table of the partitions (You need to click the **show/hide
code output** link to see it). The membership matrix (columns with name `p*`)
is inferred by
[`clue::cl_consensus()`](https://www.rdocumentation.org/link/cl_consensus?package=clue)
function with the `SE` method. Basically the value in the membership matrix
represents the probability to belong to a certain group. The finall subgroup
label for an item is determined with the group with highest probability it
belongs to.

In `get_classes()` function, the entropy is calculated from the membership
matrix and the silhouette score is calculated from the consensus matrix.



<script>
$( function() {
	$( '#tabs-node-01-get-classes' ).tabs();
} );
</script>
<div id='tabs-node-01-get-classes'>
<ul>
<li><a href='#tab-node-01-get-classes-1'>k = 2</a></li>
<li><a href='#tab-node-01-get-classes-2'>k = 3</a></li>
<li><a href='#tab-node-01-get-classes-3'>k = 4</a></li>
</ul>

<div id='tab-node-01-get-classes-1'>
<p><a id='tab-node-01-get-classes-1-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 2), get_membership(res, k = 2))
</code></pre>

<pre><code>#&gt;           class entropy silhouette   p1   p2
#&gt; sample_42     2   0.000      0.980 0.00 1.00
#&gt; sample_47     1   0.000      0.998 1.00 0.00
#&gt; sample_48     1   0.000      0.998 1.00 0.00
#&gt; sample_41     1   0.000      0.998 1.00 0.00
#&gt; sample_43     2   0.000      0.980 0.00 1.00
#&gt; sample_44     2   0.469      0.883 0.10 0.90
#&gt; sample_45     2   0.881      0.585 0.30 0.70
#&gt; sample_46     2   0.000      0.980 0.00 1.00
#&gt; sample_70     2   0.000      0.980 0.00 1.00
#&gt; sample_71     2   0.000      0.980 0.00 1.00
#&gt; sample_72     2   0.000      0.980 0.00 1.00
#&gt; sample_68     1   0.000      0.998 1.00 0.00
#&gt; sample_69     1   0.000      0.998 1.00 0.00
#&gt; sample_67     2   0.000      0.980 0.00 1.00
#&gt; sample_59     2   0.000      0.980 0.00 1.00
#&gt; sample_54     2   0.000      0.980 0.00 1.00
#&gt; sample_60     2   0.000      0.980 0.00 1.00
#&gt; sample_66     2   0.000      0.980 0.00 1.00
#&gt; sample_2      2   0.000      0.980 0.00 1.00
#&gt; sample_5      1   0.000      0.998 1.00 0.00
#&gt; sample_9      2   0.000      0.980 0.00 1.00
#&gt; sample_10     2   0.000      0.980 0.00 1.00
#&gt; sample_11     2   0.000      0.980 0.00 1.00
#&gt; sample_13     1   0.000      0.998 1.00 0.00
#&gt; sample_14     1   0.141      0.979 0.98 0.02
#&gt; sample_15     1   0.000      0.998 1.00 0.00
#&gt; sample_16     1   0.000      0.998 1.00 0.00
#&gt; sample_17     1   0.000      0.998 1.00 0.00
#&gt; sample_18     2   0.000      0.980 0.00 1.00
#&gt; sample_19     2   0.000      0.980 0.00 1.00
#&gt; sample_20     1   0.000      0.998 1.00 0.00
#&gt; sample_21     1   0.000      0.998 1.00 0.00
#&gt; sample_24     1   0.000      0.998 1.00 0.00
#&gt; sample_26     2   0.000      0.980 0.00 1.00
#&gt; sample_29     2   0.000      0.980 0.00 1.00
</code></pre>

<script>
$('#tab-node-01-get-classes-1-a').parent().next().next().hide();
$('#tab-node-01-get-classes-1-a').click(function(){
  $('#tab-node-01-get-classes-1-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-01-get-classes-2'>
<p><a id='tab-node-01-get-classes-2-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 3), get_membership(res, k = 3))
</code></pre>

<pre><code>#&gt;           class entropy silhouette   p1   p2  p3
#&gt; sample_42     3   0.000      0.990 0.00 0.00 1.0
#&gt; sample_47     2   0.628      0.159 0.46 0.54 0.0
#&gt; sample_48     1   0.000      1.000 1.00 0.00 0.0
#&gt; sample_41     1   0.000      1.000 1.00 0.00 0.0
#&gt; sample_43     2   0.000      0.958 0.00 1.00 0.0
#&gt; sample_44     2   0.000      0.958 0.00 1.00 0.0
#&gt; sample_45     2   0.000      0.958 0.00 1.00 0.0
#&gt; sample_46     2   0.000      0.958 0.00 1.00 0.0
#&gt; sample_70     2   0.000      0.958 0.00 1.00 0.0
#&gt; sample_71     3   0.000      0.990 0.00 0.00 1.0
#&gt; sample_72     3   0.000      0.990 0.00 0.00 1.0
#&gt; sample_68     1   0.000      1.000 1.00 0.00 0.0
#&gt; sample_69     1   0.000      1.000 1.00 0.00 0.0
#&gt; sample_67     3   0.000      0.990 0.00 0.00 1.0
#&gt; sample_59     2   0.000      0.958 0.00 1.00 0.0
#&gt; sample_54     2   0.000      0.958 0.00 1.00 0.0
#&gt; sample_60     2   0.000      0.958 0.00 1.00 0.0
#&gt; sample_66     3   0.000      0.990 0.00 0.00 1.0
#&gt; sample_2      3   0.000      0.990 0.00 0.00 1.0
#&gt; sample_5      1   0.000      1.000 1.00 0.00 0.0
#&gt; sample_9      3   0.000      0.990 0.00 0.00 1.0
#&gt; sample_10     3   0.000      0.990 0.00 0.00 1.0
#&gt; sample_11     3   0.000      0.990 0.00 0.00 1.0
#&gt; sample_13     1   0.000      1.000 1.00 0.00 0.0
#&gt; sample_14     3   0.296      0.888 0.10 0.00 0.9
#&gt; sample_15     1   0.000      1.000 1.00 0.00 0.0
#&gt; sample_16     2   0.153      0.924 0.04 0.96 0.0
#&gt; sample_17     1   0.000      1.000 1.00 0.00 0.0
#&gt; sample_18     2   0.000      0.958 0.00 1.00 0.0
#&gt; sample_19     2   0.000      0.958 0.00 1.00 0.0
#&gt; sample_20     1   0.000      1.000 1.00 0.00 0.0
#&gt; sample_21     1   0.000      1.000 1.00 0.00 0.0
#&gt; sample_24     1   0.000      1.000 1.00 0.00 0.0
#&gt; sample_26     2   0.000      0.958 0.00 1.00 0.0
#&gt; sample_29     3   0.000      0.990 0.00 0.00 1.0
</code></pre>

<script>
$('#tab-node-01-get-classes-2-a').parent().next().next().hide();
$('#tab-node-01-get-classes-2-a').click(function(){
  $('#tab-node-01-get-classes-2-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-01-get-classes-3'>
<p><a id='tab-node-01-get-classes-3-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 4), get_membership(res, k = 4))
</code></pre>

<pre><code>#&gt;           class entropy silhouette   p1   p2   p3   p4
#&gt; sample_42     3   0.000      0.921 0.00 0.00 1.00 0.00
#&gt; sample_47     4   0.648      0.426 0.30 0.10 0.00 0.60
#&gt; sample_48     1   0.000      1.000 1.00 0.00 0.00 0.00
#&gt; sample_41     1   0.000      1.000 1.00 0.00 0.00 0.00
#&gt; sample_43     2   0.000      0.801 0.00 1.00 0.00 0.00
#&gt; sample_44     2   0.000      0.801 0.00 1.00 0.00 0.00
#&gt; sample_45     2   0.000      0.801 0.00 1.00 0.00 0.00
#&gt; sample_46     2   0.000      0.801 0.00 1.00 0.00 0.00
#&gt; sample_70     2   0.000      0.801 0.00 1.00 0.00 0.00
#&gt; sample_71     3   0.000      0.921 0.00 0.00 1.00 0.00
#&gt; sample_72     3   0.000      0.921 0.00 0.00 1.00 0.00
#&gt; sample_68     1   0.000      1.000 1.00 0.00 0.00 0.00
#&gt; sample_69     1   0.000      1.000 1.00 0.00 0.00 0.00
#&gt; sample_67     3   0.000      0.921 0.00 0.00 1.00 0.00
#&gt; sample_59     2   0.000      0.801 0.00 1.00 0.00 0.00
#&gt; sample_54     2   0.499      0.226 0.00 0.52 0.00 0.48
#&gt; sample_60     2   0.499      0.226 0.00 0.52 0.00 0.48
#&gt; sample_66     3   0.265      0.893 0.00 0.00 0.88 0.12
#&gt; sample_2      3   0.121      0.916 0.00 0.00 0.96 0.04
#&gt; sample_5      1   0.000      1.000 1.00 0.00 0.00 0.00
#&gt; sample_9      3   0.441      0.714 0.00 0.00 0.70 0.30
#&gt; sample_10     3   0.317      0.867 0.00 0.00 0.84 0.16
#&gt; sample_11     3   0.265      0.893 0.00 0.00 0.88 0.12
#&gt; sample_13     1   0.000      1.000 1.00 0.00 0.00 0.00
#&gt; sample_14     4   0.543     -0.119 0.02 0.00 0.38 0.60
#&gt; sample_15     1   0.000      1.000 1.00 0.00 0.00 0.00
#&gt; sample_16     4   0.380      0.435 0.00 0.22 0.00 0.78
#&gt; sample_17     4   0.234      0.542 0.10 0.00 0.00 0.90
#&gt; sample_18     2   0.471      0.461 0.00 0.64 0.00 0.36
#&gt; sample_19     4   0.462      0.225 0.00 0.34 0.00 0.66
#&gt; sample_20     1   0.000      1.000 1.00 0.00 0.00 0.00
#&gt; sample_21     1   0.000      1.000 1.00 0.00 0.00 0.00
#&gt; sample_24     1   0.000      1.000 1.00 0.00 0.00 0.00
#&gt; sample_26     2   0.234      0.745 0.00 0.90 0.00 0.10
#&gt; sample_29     3   0.000      0.921 0.00 0.00 1.00 0.00
</code></pre>

<script>
$('#tab-node-01-get-classes-3-a').parent().next().next().hide();
$('#tab-node-01-get-classes-3-a').click(function(){
  $('#tab-node-01-get-classes-3-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>
</div>

Heatmaps for the consensus matrix. It visualizes the probability of two
samples to be in a same group.




<script>
$( function() {
	$( '#tabs-node-01-consensus-heatmap' ).tabs();
} );
</script>
<div id='tabs-node-01-consensus-heatmap'>
<ul>
<li><a href='#tab-node-01-consensus-heatmap-1'>k = 2</a></li>
<li><a href='#tab-node-01-consensus-heatmap-2'>k = 3</a></li>
<li><a href='#tab-node-01-consensus-heatmap-3'>k = 4</a></li>
</ul>
<div id='tab-node-01-consensus-heatmap-1'>
<pre><code class="r">consensus_heatmap(res, k = 2)
</code></pre>

<p><img src="figure_cola/tab-node-01-consensus-heatmap-1-1.png" alt="plot of chunk tab-node-01-consensus-heatmap-1"/></p>

</div>
<div id='tab-node-01-consensus-heatmap-2'>
<pre><code class="r">consensus_heatmap(res, k = 3)
</code></pre>

<p><img src="figure_cola/tab-node-01-consensus-heatmap-2-1.png" alt="plot of chunk tab-node-01-consensus-heatmap-2"/></p>

</div>
<div id='tab-node-01-consensus-heatmap-3'>
<pre><code class="r">consensus_heatmap(res, k = 4)
</code></pre>

<p><img src="figure_cola/tab-node-01-consensus-heatmap-3-1.png" alt="plot of chunk tab-node-01-consensus-heatmap-3"/></p>

</div>
</div>

Heatmaps for the membership of samples in all partitions to see how consistent they are:





<script>
$( function() {
	$( '#tabs-node-01-membership-heatmap' ).tabs();
} );
</script>
<div id='tabs-node-01-membership-heatmap'>
<ul>
<li><a href='#tab-node-01-membership-heatmap-1'>k = 2</a></li>
<li><a href='#tab-node-01-membership-heatmap-2'>k = 3</a></li>
<li><a href='#tab-node-01-membership-heatmap-3'>k = 4</a></li>
</ul>
<div id='tab-node-01-membership-heatmap-1'>
<pre><code class="r">membership_heatmap(res, k = 2)
</code></pre>

<p><img src="figure_cola/tab-node-01-membership-heatmap-1-1.png" alt="plot of chunk tab-node-01-membership-heatmap-1"/></p>

</div>
<div id='tab-node-01-membership-heatmap-2'>
<pre><code class="r">membership_heatmap(res, k = 3)
</code></pre>

<p><img src="figure_cola/tab-node-01-membership-heatmap-2-1.png" alt="plot of chunk tab-node-01-membership-heatmap-2"/></p>

</div>
<div id='tab-node-01-membership-heatmap-3'>
<pre><code class="r">membership_heatmap(res, k = 4)
</code></pre>

<p><img src="figure_cola/tab-node-01-membership-heatmap-3-1.png" alt="plot of chunk tab-node-01-membership-heatmap-3"/></p>

</div>
</div>

As soon as the classes for columns are determined, the signatures
that are significantly different between subgroups can be looked for. 
Following are the heatmaps for signatures.




Signature heatmaps where rows are scaled:



<script>
$( function() {
	$( '#tabs-node-01-get-signatures' ).tabs();
} );
</script>
<div id='tabs-node-01-get-signatures'>
<ul>
<li><a href='#tab-node-01-get-signatures-1'>k = 2</a></li>
<li><a href='#tab-node-01-get-signatures-2'>k = 3</a></li>
<li><a href='#tab-node-01-get-signatures-3'>k = 4</a></li>
</ul>
<div id='tab-node-01-get-signatures-1'>
<pre><code class="r">get_signatures(res, k = 2)
</code></pre>

<p><img src="figure_cola/tab-node-01-get-signatures-1-1.png" alt="plot of chunk tab-node-01-get-signatures-1"/></p>

</div>
<div id='tab-node-01-get-signatures-2'>
<pre><code class="r">get_signatures(res, k = 3)
</code></pre>

<p><img src="figure_cola/tab-node-01-get-signatures-2-1.png" alt="plot of chunk tab-node-01-get-signatures-2"/></p>

</div>
<div id='tab-node-01-get-signatures-3'>
<pre><code class="r">get_signatures(res, k = 4)
</code></pre>

<p><img src="figure_cola/tab-node-01-get-signatures-3-1.png" alt="plot of chunk tab-node-01-get-signatures-3"/></p>

</div>
</div>



Signature heatmaps where rows are not scaled:


<script>
$( function() {
	$( '#tabs-node-01-get-signatures-no-scale' ).tabs();
} );
</script>
<div id='tabs-node-01-get-signatures-no-scale'>
<ul>
<li><a href='#tab-node-01-get-signatures-no-scale-1'>k = 2</a></li>
<li><a href='#tab-node-01-get-signatures-no-scale-2'>k = 3</a></li>
<li><a href='#tab-node-01-get-signatures-no-scale-3'>k = 4</a></li>
</ul>
<div id='tab-node-01-get-signatures-no-scale-1'>
<pre><code class="r">get_signatures(res, k = 2, scale_rows = FALSE)
</code></pre>

<p><img src="figure_cola/tab-node-01-get-signatures-no-scale-1-1.png" alt="plot of chunk tab-node-01-get-signatures-no-scale-1"/></p>

</div>
<div id='tab-node-01-get-signatures-no-scale-2'>
<pre><code class="r">get_signatures(res, k = 3, scale_rows = FALSE)
</code></pre>

<p><img src="figure_cola/tab-node-01-get-signatures-no-scale-2-1.png" alt="plot of chunk tab-node-01-get-signatures-no-scale-2"/></p>

</div>
<div id='tab-node-01-get-signatures-no-scale-3'>
<pre><code class="r">get_signatures(res, k = 4, scale_rows = FALSE)
</code></pre>

<p><img src="figure_cola/tab-node-01-get-signatures-no-scale-3-1.png" alt="plot of chunk tab-node-01-get-signatures-no-scale-3"/></p>

</div>
</div>



Compare the overlap of signatures from different k:

```r
compare_signatures(res)
```

![plot of chunk node-01-signature_compare](figure_cola/node-01-signature_compare-1.png)

`get_signature()` returns a data frame invisibly. To get the list of signatures, the function
call should be assigned to a variable explicitly. In following code, if `plot` argument is set
to `FALSE`, no heatmap is plotted while only the differential analysis is performed.

```r
# code only for demonstration
tb = get_signature(res, k = ..., plot = FALSE)
```

An example of the output of `tb` is:

```
#>   which_row         fdr    mean_1    mean_2 scaled_mean_1 scaled_mean_2 km
#> 1        38 0.042760348  8.373488  9.131774    -0.5533452     0.5164555  1
#> 2        40 0.018707592  7.106213  8.469186    -0.6173731     0.5762149  1
#> 3        55 0.019134737 10.221463 11.207825    -0.6159697     0.5749050  1
#> 4        59 0.006059896  5.921854  7.869574    -0.6899429     0.6439467  1
#> 5        60 0.018055526  8.928898 10.211722    -0.6204761     0.5791110  1
#> 6        98 0.009384629 15.714769 14.887706     0.6635654    -0.6193277  2
...
```

The columns in `tb` are:

1. `which_row`: row indices corresponding to the input matrix.
2. `fdr`: FDR for the differential test. 
3. `mean_x`: The mean value in group x.
4. `scaled_mean_x`: The mean value in group x after rows are scaled.
5. `km`: Row groups if k-means clustering is applied to rows (which is done by automatically selecting number of clusters).

If there are too many signatures, `top_signatures = ...` can be set to only show the 
signatures with the highest FDRs:

```r
# code only for demonstration
# e.g. to show the top 500 most significant rows
tb = get_signature(res, k = ..., top_signatures = 500)
```

If the signatures are defined as these which are uniquely high in current group, `diff_method` argument
can be set to `"uniquely_high_in_one_group"`:

```r
# code only for demonstration
tb = get_signature(res, k = ..., diff_method = "uniquely_high_in_one_group")
```




UMAP plot which shows how samples are separated.


<script>
$( function() {
	$( '#tabs-node-01-dimension-reduction' ).tabs();
} );
</script>
<div id='tabs-node-01-dimension-reduction'>
<ul>
<li><a href='#tab-node-01-dimension-reduction-1'>k = 2</a></li>
<li><a href='#tab-node-01-dimension-reduction-2'>k = 3</a></li>
<li><a href='#tab-node-01-dimension-reduction-3'>k = 4</a></li>
</ul>
<div id='tab-node-01-dimension-reduction-1'>
<pre><code class="r">dimension_reduction(res, k = 2, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-01-dimension-reduction-1-1.png" alt="plot of chunk tab-node-01-dimension-reduction-1"/></p>

</div>
<div id='tab-node-01-dimension-reduction-2'>
<pre><code class="r">dimension_reduction(res, k = 3, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-01-dimension-reduction-2-1.png" alt="plot of chunk tab-node-01-dimension-reduction-2"/></p>

</div>
<div id='tab-node-01-dimension-reduction-3'>
<pre><code class="r">dimension_reduction(res, k = 4, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-01-dimension-reduction-3-1.png" alt="plot of chunk tab-node-01-dimension-reduction-3"/></p>

</div>
</div>



Following heatmap shows how subgroups are split when increasing `k`:

```r
collect_classes(res)
```

![plot of chunk node-01-collect-classes](figure_cola/node-01-collect-classes-1.png)




Test correlation between subgroups and known annotations. If the known
annotation is numeric, one-way ANOVA test is applied, and if the known
annotation is discrete, chi-squared contingency table test is applied.

```r
test_to_known_factors(res)
```

```
#>             n_sample ALL.AML(p-value) k
#> ATC:skmeans       35            0.233 2
#> ATC:skmeans       34            0.336 3
#> ATC:skmeans       28            0.275 4
```




If matrix rows can be associated to genes, consider to use `functional_enrichment(res,
...)` to perform function enrichment for the signature genes. See [this vignette](https://jokergoo.github.io/cola_vignettes/functional_enrichment.html) for more detailed explanations.


 

---------------------------------------------------




### Node02


Parent node: [Node0](#Node0).
Child nodes: 
                Node011-leaf
        ,
                Node012-leaf
        ,
                Node013-leaf
        ,
                Node021-leaf
        ,
                Node022-leaf
        .







The object with results only for a single top-value method and a single partitioning method 
can be extracted as:

```r
res = res_rh["02"]
```

A summary of `res` and all the functions that can be applied to it:

```r
res
```

```
#> A 'ConsensusPartition' object with k = 2, 3, 4.
#>   On a matrix with 3910 rows and 24 columns.
#>   Top rows (391) are extracted by 'ATC' method.
#>   Subgroups are detected by 'skmeans' method.
#>   Performed in total 150 partitions by row resampling.
#>   Best k for subgroups seems to be 2.
#> 
#> Following methods can be applied to this 'ConsensusPartition' object:
#>  [1] "cola_report"             "collect_classes"         "collect_plots"          
#>  [4] "collect_stats"           "colnames"                "compare_partitions"     
#>  [7] "compare_signatures"      "consensus_heatmap"       "dimension_reduction"    
#> [10] "functional_enrichment"   "get_anno_col"            "get_anno"               
#> [13] "get_classes"             "get_consensus"           "get_matrix"             
#> [16] "get_membership"          "get_param"               "get_signatures"         
#> [19] "get_stats"               "is_best_k"               "is_stable_k"            
#> [22] "membership_heatmap"      "ncol"                    "nrow"                   
#> [25] "plot_ecdf"               "predict_classes"         "rownames"               
#> [28] "select_partition_number" "show"                    "suggest_best_k"         
#> [31] "test_to_known_factors"   "top_rows_heatmap"
```

`collect_plots()` function collects all the plots made from `res` for all `k` (number of subgroups)
into one single page to provide an easy and fast comparison between different `k`.

```r
collect_plots(res)
```

![plot of chunk node-02-collect-plots](figure_cola/node-02-collect-plots-1.png)

The plots are:

- The first row: a plot of the eCDF (empirical cumulative distribution
  function) curves of the consensus matrix for each `k` and the heatmap of
  predicted classes for each `k`.
- The second row: heatmaps of the consensus matrix for each `k`.
- The third row: heatmaps of the membership matrix for each `k`.
- The fouth row: heatmaps of the signatures for each `k`.

All the plots in panels can be made by individual functions and they are
plotted later in this section.

`select_partition_number()` produces several plots showing different
statistics for choosing "optimized" `k`. There are following statistics:

- eCDF curves of the consensus matrix for each `k`;
- 1-PAC. [The PAC score](https://en.wikipedia.org/wiki/Consensus_clustering#Over-interpretation_potential_of_consensus_clustering)
  measures the proportion of the ambiguous subgrouping.
- Mean silhouette score.
- Concordance. The mean probability of fiting the consensus subgroup labels in all
  partitions.
- Area increased. Denote $A_k$ as the area under the eCDF curve for current
  `k`, the area increased is defined as $A_k - A_{k-1}$.
- Rand index. The percent of pairs of samples that are both in a same cluster
  or both are not in a same cluster in the partition of k and k-1.
- Jaccard index. The ratio of pairs of samples are both in a same cluster in
  the partition of k and k-1 and the pairs of samples are both in a same
  cluster in the partition k or k-1.

The detailed explanations of these statistics can be found in [the _cola_
vignette](https://jokergoo.github.io/cola_vignettes/cola.html#toc_13).

Generally speaking, higher 1-PAC score, higher mean silhouette score or higher
concordance corresponds to better partition. Rand index and Jaccard index
measure how similar the current partition is compared to partition with `k-1`.
If they are too similar, we won't accept `k` is better than `k-1`.

```r
select_partition_number(res)
```

![plot of chunk node-02-select-partition-number](figure_cola/node-02-select-partition-number-1.png)

The numeric values for all these statistics can be obtained by `get_stats()`.

```r
get_stats(res)
```

```
#>   k 1-PAC mean_silhouette concordance area_increased  Rand Jaccard
#> 2 2 1.000           1.000       1.000          0.519 0.482   0.482
#> 3 3 0.823           0.949       0.966          0.294 0.819   0.635
#> 4 4 0.870           0.880       0.928          0.128 0.909   0.725
```

`suggest_best_k()` suggests the best $k$ based on these statistics. The rules are as follows:

- All $k$ with Jaccard index larger than 0.95 are removed because increasing
  $k$ does not provide enough extra information. If all $k$ are removed, it is
  marked as no subgroup is detected.
- For all $k$ with 1-PAC score larger than 0.9, the maximal $k$ is taken as
  the best $k$, and other $k$ are marked as optional $k$.
- If it does not fit the second rule. The $k$ with the maximal vote of the
  highest 1-PAC score, highest mean silhouette, and highest concordance is
  taken as the best $k$.

```r
suggest_best_k(res)
```

```
#> [1] 2
```


Following is the table of the partitions (You need to click the **show/hide
code output** link to see it). The membership matrix (columns with name `p*`)
is inferred by
[`clue::cl_consensus()`](https://www.rdocumentation.org/link/cl_consensus?package=clue)
function with the `SE` method. Basically the value in the membership matrix
represents the probability to belong to a certain group. The finall subgroup
label for an item is determined with the group with highest probability it
belongs to.

In `get_classes()` function, the entropy is calculated from the membership
matrix and the silhouette score is calculated from the consensus matrix.



<script>
$( function() {
	$( '#tabs-node-02-get-classes' ).tabs();
} );
</script>
<div id='tabs-node-02-get-classes'>
<ul>
<li><a href='#tab-node-02-get-classes-1'>k = 2</a></li>
<li><a href='#tab-node-02-get-classes-2'>k = 3</a></li>
<li><a href='#tab-node-02-get-classes-3'>k = 4</a></li>
</ul>

<div id='tab-node-02-get-classes-1'>
<p><a id='tab-node-02-get-classes-1-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 2), get_membership(res, k = 2))
</code></pre>

<pre><code>#&gt;           class entropy silhouette p1 p2
#&gt; sample_52     2       0          1  0  1
#&gt; sample_53     1       0          1  1  0
#&gt; sample_51     1       0          1  1  0
#&gt; sample_50     1       0          1  1  0
#&gt; sample_57     2       0          1  0  1
#&gt; sample_58     2       0          1  0  1
#&gt; sample_61     2       0          1  0  1
#&gt; sample_65     2       0          1  0  1
#&gt; sample_63     2       0          1  0  1
#&gt; sample_64     2       0          1  0  1
#&gt; sample_62     2       0          1  0  1
#&gt; sample_12     1       0          1  1  0
#&gt; sample_22     1       0          1  1  0
#&gt; sample_25     1       0          1  1  0
#&gt; sample_34     2       0          1  0  1
#&gt; sample_35     2       0          1  0  1
#&gt; sample_36     1       0          1  1  0
#&gt; sample_37     1       0          1  1  0
#&gt; sample_38     1       0          1  1  0
#&gt; sample_28     1       0          1  1  0
#&gt; sample_30     1       0          1  1  0
#&gt; sample_31     2       0          1  0  1
#&gt; sample_32     1       0          1  1  0
#&gt; sample_33     1       0          1  1  0
</code></pre>

<script>
$('#tab-node-02-get-classes-1-a').parent().next().next().hide();
$('#tab-node-02-get-classes-1-a').click(function(){
  $('#tab-node-02-get-classes-1-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-02-get-classes-2'>
<p><a id='tab-node-02-get-classes-2-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 3), get_membership(res, k = 3))
</code></pre>

<pre><code>#&gt;           class entropy silhouette   p1   p2   p3
#&gt; sample_52     2  0.2959      0.931 0.00 0.90 0.10
#&gt; sample_53     1  0.0892      0.981 0.98 0.00 0.02
#&gt; sample_51     1  0.0892      0.981 0.98 0.00 0.02
#&gt; sample_50     3  0.2959      0.891 0.10 0.00 0.90
#&gt; sample_57     2  0.0000      0.956 0.00 1.00 0.00
#&gt; sample_58     2  0.0000      0.956 0.00 1.00 0.00
#&gt; sample_61     2  0.0000      0.956 0.00 1.00 0.00
#&gt; sample_65     2  0.2959      0.931 0.00 0.90 0.10
#&gt; sample_63     2  0.2959      0.931 0.00 0.90 0.10
#&gt; sample_64     2  0.0000      0.956 0.00 1.00 0.00
#&gt; sample_62     2  0.2959      0.931 0.00 0.90 0.10
#&gt; sample_12     1  0.2537      0.907 0.92 0.00 0.08
#&gt; sample_22     1  0.0000      0.976 1.00 0.00 0.00
#&gt; sample_25     1  0.0000      0.976 1.00 0.00 0.00
#&gt; sample_34     2  0.0000      0.956 0.00 1.00 0.00
#&gt; sample_35     2  0.0000      0.956 0.00 1.00 0.00
#&gt; sample_36     3  0.0000      0.933 0.00 0.00 1.00
#&gt; sample_37     1  0.0892      0.981 0.98 0.00 0.02
#&gt; sample_38     1  0.0000      0.976 1.00 0.00 0.00
#&gt; sample_28     1  0.0892      0.981 0.98 0.00 0.02
#&gt; sample_30     3  0.0000      0.933 0.00 0.00 1.00
#&gt; sample_31     3  0.0892      0.919 0.00 0.02 0.98
#&gt; sample_32     3  0.3832      0.888 0.10 0.02 0.88
#&gt; sample_33     1  0.0892      0.981 0.98 0.00 0.02
</code></pre>

<script>
$('#tab-node-02-get-classes-2-a').parent().next().next().hide();
$('#tab-node-02-get-classes-2-a').click(function(){
  $('#tab-node-02-get-classes-2-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>

<div id='tab-node-02-get-classes-3'>
<p><a id='tab-node-02-get-classes-3-a' style='color:#0366d6' href='#'>show/hide code output</a></p>
<pre><code class="r">cbind(get_classes(res, k = 4), get_membership(res, k = 4))
</code></pre>

<pre><code>#&gt;           class entropy silhouette   p1   p2   p3   p4
#&gt; sample_52     2   0.000      0.987 0.00 1.00 0.00 0.00
#&gt; sample_53     1   0.000      0.957 1.00 0.00 0.00 0.00
#&gt; sample_51     1   0.000      0.957 1.00 0.00 0.00 0.00
#&gt; sample_50     3   0.000      0.900 0.00 0.00 1.00 0.00
#&gt; sample_57     2   0.121      0.945 0.00 0.96 0.00 0.04
#&gt; sample_58     4   0.499      0.277 0.00 0.48 0.00 0.52
#&gt; sample_61     4   0.201      0.856 0.00 0.08 0.00 0.92
#&gt; sample_65     2   0.000      0.987 0.00 1.00 0.00 0.00
#&gt; sample_63     2   0.000      0.987 0.00 1.00 0.00 0.00
#&gt; sample_64     4   0.292      0.848 0.00 0.14 0.00 0.86
#&gt; sample_62     2   0.000      0.987 0.00 1.00 0.00 0.00
#&gt; sample_12     1   0.484      0.795 0.78 0.00 0.14 0.08
#&gt; sample_22     1   0.201      0.926 0.92 0.00 0.00 0.08
#&gt; sample_25     1   0.234      0.916 0.90 0.00 0.00 0.10
#&gt; sample_34     4   0.201      0.857 0.00 0.08 0.00 0.92
#&gt; sample_35     4   0.234      0.863 0.00 0.10 0.00 0.90
#&gt; sample_36     3   0.000      0.900 0.00 0.00 1.00 0.00
#&gt; sample_37     1   0.000      0.957 1.00 0.00 0.00 0.00
#&gt; sample_38     1   0.000      0.957 1.00 0.00 0.00 0.00
#&gt; sample_28     1   0.000      0.957 1.00 0.00 0.00 0.00
#&gt; sample_30     3   0.000      0.900 0.00 0.00 1.00 0.00
#&gt; sample_31     3   0.506      0.570 0.00 0.30 0.68 0.02
#&gt; sample_32     3   0.261      0.867 0.02 0.02 0.92 0.04
#&gt; sample_33     1   0.000      0.957 1.00 0.00 0.00 0.00
</code></pre>

<script>
$('#tab-node-02-get-classes-3-a').parent().next().next().hide();
$('#tab-node-02-get-classes-3-a').click(function(){
  $('#tab-node-02-get-classes-3-a').parent().next().next().toggle();
  return(false);
});
</script>
</div>
</div>

Heatmaps for the consensus matrix. It visualizes the probability of two
samples to be in a same group.




<script>
$( function() {
	$( '#tabs-node-02-consensus-heatmap' ).tabs();
} );
</script>
<div id='tabs-node-02-consensus-heatmap'>
<ul>
<li><a href='#tab-node-02-consensus-heatmap-1'>k = 2</a></li>
<li><a href='#tab-node-02-consensus-heatmap-2'>k = 3</a></li>
<li><a href='#tab-node-02-consensus-heatmap-3'>k = 4</a></li>
</ul>
<div id='tab-node-02-consensus-heatmap-1'>
<pre><code class="r">consensus_heatmap(res, k = 2)
</code></pre>

<p><img src="figure_cola/tab-node-02-consensus-heatmap-1-1.png" alt="plot of chunk tab-node-02-consensus-heatmap-1"/></p>

</div>
<div id='tab-node-02-consensus-heatmap-2'>
<pre><code class="r">consensus_heatmap(res, k = 3)
</code></pre>

<p><img src="figure_cola/tab-node-02-consensus-heatmap-2-1.png" alt="plot of chunk tab-node-02-consensus-heatmap-2"/></p>

</div>
<div id='tab-node-02-consensus-heatmap-3'>
<pre><code class="r">consensus_heatmap(res, k = 4)
</code></pre>

<p><img src="figure_cola/tab-node-02-consensus-heatmap-3-1.png" alt="plot of chunk tab-node-02-consensus-heatmap-3"/></p>

</div>
</div>

Heatmaps for the membership of samples in all partitions to see how consistent they are:





<script>
$( function() {
	$( '#tabs-node-02-membership-heatmap' ).tabs();
} );
</script>
<div id='tabs-node-02-membership-heatmap'>
<ul>
<li><a href='#tab-node-02-membership-heatmap-1'>k = 2</a></li>
<li><a href='#tab-node-02-membership-heatmap-2'>k = 3</a></li>
<li><a href='#tab-node-02-membership-heatmap-3'>k = 4</a></li>
</ul>
<div id='tab-node-02-membership-heatmap-1'>
<pre><code class="r">membership_heatmap(res, k = 2)
</code></pre>

<p><img src="figure_cola/tab-node-02-membership-heatmap-1-1.png" alt="plot of chunk tab-node-02-membership-heatmap-1"/></p>

</div>
<div id='tab-node-02-membership-heatmap-2'>
<pre><code class="r">membership_heatmap(res, k = 3)
</code></pre>

<p><img src="figure_cola/tab-node-02-membership-heatmap-2-1.png" alt="plot of chunk tab-node-02-membership-heatmap-2"/></p>

</div>
<div id='tab-node-02-membership-heatmap-3'>
<pre><code class="r">membership_heatmap(res, k = 4)
</code></pre>

<p><img src="figure_cola/tab-node-02-membership-heatmap-3-1.png" alt="plot of chunk tab-node-02-membership-heatmap-3"/></p>

</div>
</div>

As soon as the classes for columns are determined, the signatures
that are significantly different between subgroups can be looked for. 
Following are the heatmaps for signatures.




Signature heatmaps where rows are scaled:



<script>
$( function() {
	$( '#tabs-node-02-get-signatures' ).tabs();
} );
</script>
<div id='tabs-node-02-get-signatures'>
<ul>
<li><a href='#tab-node-02-get-signatures-1'>k = 2</a></li>
<li><a href='#tab-node-02-get-signatures-2'>k = 3</a></li>
<li><a href='#tab-node-02-get-signatures-3'>k = 4</a></li>
</ul>
<div id='tab-node-02-get-signatures-1'>
<pre><code class="r">get_signatures(res, k = 2)
</code></pre>

<p><img src="figure_cola/tab-node-02-get-signatures-1-1.png" alt="plot of chunk tab-node-02-get-signatures-1"/></p>

</div>
<div id='tab-node-02-get-signatures-2'>
<pre><code class="r">get_signatures(res, k = 3)
</code></pre>

<p><img src="figure_cola/tab-node-02-get-signatures-2-1.png" alt="plot of chunk tab-node-02-get-signatures-2"/></p>

</div>
<div id='tab-node-02-get-signatures-3'>
<pre><code class="r">get_signatures(res, k = 4)
</code></pre>

<p><img src="figure_cola/tab-node-02-get-signatures-3-1.png" alt="plot of chunk tab-node-02-get-signatures-3"/></p>

</div>
</div>



Signature heatmaps where rows are not scaled:


<script>
$( function() {
	$( '#tabs-node-02-get-signatures-no-scale' ).tabs();
} );
</script>
<div id='tabs-node-02-get-signatures-no-scale'>
<ul>
<li><a href='#tab-node-02-get-signatures-no-scale-1'>k = 2</a></li>
<li><a href='#tab-node-02-get-signatures-no-scale-2'>k = 3</a></li>
<li><a href='#tab-node-02-get-signatures-no-scale-3'>k = 4</a></li>
</ul>
<div id='tab-node-02-get-signatures-no-scale-1'>
<pre><code class="r">get_signatures(res, k = 2, scale_rows = FALSE)
</code></pre>

<p><img src="figure_cola/tab-node-02-get-signatures-no-scale-1-1.png" alt="plot of chunk tab-node-02-get-signatures-no-scale-1"/></p>

</div>
<div id='tab-node-02-get-signatures-no-scale-2'>
<pre><code class="r">get_signatures(res, k = 3, scale_rows = FALSE)
</code></pre>

<p><img src="figure_cola/tab-node-02-get-signatures-no-scale-2-1.png" alt="plot of chunk tab-node-02-get-signatures-no-scale-2"/></p>

</div>
<div id='tab-node-02-get-signatures-no-scale-3'>
<pre><code class="r">get_signatures(res, k = 4, scale_rows = FALSE)
</code></pre>

<p><img src="figure_cola/tab-node-02-get-signatures-no-scale-3-1.png" alt="plot of chunk tab-node-02-get-signatures-no-scale-3"/></p>

</div>
</div>



Compare the overlap of signatures from different k:

```r
compare_signatures(res)
```

![plot of chunk node-02-signature_compare](figure_cola/node-02-signature_compare-1.png)

`get_signature()` returns a data frame invisibly. To get the list of signatures, the function
call should be assigned to a variable explicitly. In following code, if `plot` argument is set
to `FALSE`, no heatmap is plotted while only the differential analysis is performed.

```r
# code only for demonstration
tb = get_signature(res, k = ..., plot = FALSE)
```

An example of the output of `tb` is:

```
#>   which_row         fdr    mean_1    mean_2 scaled_mean_1 scaled_mean_2 km
#> 1        38 0.042760348  8.373488  9.131774    -0.5533452     0.5164555  1
#> 2        40 0.018707592  7.106213  8.469186    -0.6173731     0.5762149  1
#> 3        55 0.019134737 10.221463 11.207825    -0.6159697     0.5749050  1
#> 4        59 0.006059896  5.921854  7.869574    -0.6899429     0.6439467  1
#> 5        60 0.018055526  8.928898 10.211722    -0.6204761     0.5791110  1
#> 6        98 0.009384629 15.714769 14.887706     0.6635654    -0.6193277  2
...
```

The columns in `tb` are:

1. `which_row`: row indices corresponding to the input matrix.
2. `fdr`: FDR for the differential test. 
3. `mean_x`: The mean value in group x.
4. `scaled_mean_x`: The mean value in group x after rows are scaled.
5. `km`: Row groups if k-means clustering is applied to rows (which is done by automatically selecting number of clusters).

If there are too many signatures, `top_signatures = ...` can be set to only show the 
signatures with the highest FDRs:

```r
# code only for demonstration
# e.g. to show the top 500 most significant rows
tb = get_signature(res, k = ..., top_signatures = 500)
```

If the signatures are defined as these which are uniquely high in current group, `diff_method` argument
can be set to `"uniquely_high_in_one_group"`:

```r
# code only for demonstration
tb = get_signature(res, k = ..., diff_method = "uniquely_high_in_one_group")
```




UMAP plot which shows how samples are separated.


<script>
$( function() {
	$( '#tabs-node-02-dimension-reduction' ).tabs();
} );
</script>
<div id='tabs-node-02-dimension-reduction'>
<ul>
<li><a href='#tab-node-02-dimension-reduction-1'>k = 2</a></li>
<li><a href='#tab-node-02-dimension-reduction-2'>k = 3</a></li>
<li><a href='#tab-node-02-dimension-reduction-3'>k = 4</a></li>
</ul>
<div id='tab-node-02-dimension-reduction-1'>
<pre><code class="r">dimension_reduction(res, k = 2, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-02-dimension-reduction-1-1.png" alt="plot of chunk tab-node-02-dimension-reduction-1"/></p>

</div>
<div id='tab-node-02-dimension-reduction-2'>
<pre><code class="r">dimension_reduction(res, k = 3, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-02-dimension-reduction-2-1.png" alt="plot of chunk tab-node-02-dimension-reduction-2"/></p>

</div>
<div id='tab-node-02-dimension-reduction-3'>
<pre><code class="r">dimension_reduction(res, k = 4, method = &quot;UMAP&quot;)
</code></pre>

<p><img src="figure_cola/tab-node-02-dimension-reduction-3-1.png" alt="plot of chunk tab-node-02-dimension-reduction-3"/></p>

</div>
</div>



Following heatmap shows how subgroups are split when increasing `k`:

```r
collect_classes(res)
```

![plot of chunk node-02-collect-classes](figure_cola/node-02-collect-classes-1.png)




Test correlation between subgroups and known annotations. If the known
annotation is numeric, one-way ANOVA test is applied, and if the known
annotation is discrete, chi-squared contingency table test is applied.

```r
test_to_known_factors(res)
```

```
#>             n_sample ALL.AML(p-value) k
#> ATC:skmeans       24           0.2784 2
#> ATC:skmeans       24           0.0574 3
#> ATC:skmeans       23           0.1468 4
```




If matrix rows can be associated to genes, consider to use `functional_enrichment(res,
...)` to perform function enrichment for the signature genes. See [this vignette](https://jokergoo.github.io/cola_vignettes/functional_enrichment.html) for more detailed explanations.


 

## Session info


```r
sessionInfo()
```

```
#> R version 4.1.0 (2021-05-18)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: CentOS Linux 7 (Core)
#> 
#> Matrix products: default
#> BLAS/LAPACK: /usr/lib64/libopenblas-r0.3.3.so
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
#>  [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
#> [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] genefilter_1.74.0     ComplexHeatmap_2.8.0  markdown_1.1          knitr_1.33           
#> [5] preprocessCore_1.54.0 golubEsets_1.34.0     Biobase_2.52.0        BiocGenerics_0.38.0  
#> [9] cola_1.9.4           
#> 
#> loaded via a namespace (and not attached):
#>   [1] colorspace_2.0-2       rjson_0.2.20           ellipsis_0.3.2         mclust_5.4.7          
#>   [5] circlize_0.4.13        XVector_0.32.0         GlobalOptions_0.1.2    clue_0.3-59           
#>   [9] rstudioapi_0.13        bit64_4.0.5            AnnotationDbi_1.54.1   Polychrome_1.3.1      
#>  [13] RSpectra_0.16-0        fansi_0.5.0            xml2_1.3.2             codetools_0.2-18      
#>  [17] splines_4.1.0          doParallel_1.0.16      cachem_1.0.5           impute_1.66.0         
#>  [21] polyclip_1.10-0        jsonlite_1.7.2         Cairo_1.5-12.2         umap_0.2.7.0          
#>  [25] annotate_1.70.0        cluster_2.1.2          png_0.1-7              data.tree_1.0.0       
#>  [29] compiler_4.1.0         httr_1.4.2             assertthat_0.2.1       Matrix_1.3-4          
#>  [33] fastmap_1.1.0          tools_4.1.0            gtable_0.3.0           glue_1.4.2            
#>  [37] GenomeInfoDbData_1.2.6 dplyr_1.0.7            Rcpp_1.0.7             slam_0.1-48           
#>  [41] eulerr_6.1.0           vctrs_0.3.8            Biostrings_2.60.1      iterators_1.0.13      
#>  [45] polylabelr_0.2.0       xfun_0.24              stringr_1.4.0          lifecycle_1.0.0       
#>  [49] irlba_2.3.3            XML_3.99-0.6           dendextend_1.15.1      zlibbioc_1.38.0       
#>  [53] scales_1.1.1           microbenchmark_1.4-7   RColorBrewer_1.1-2     memoise_2.0.0         
#>  [57] reticulate_1.20        gridExtra_2.3          ggplot2_3.3.5          stringi_1.7.3         
#>  [61] RSQLite_2.2.7          highr_0.9              S4Vectors_0.30.0       foreach_1.5.1         
#>  [65] shape_1.4.6            GenomeInfoDb_1.28.1    rlang_0.4.11           pkgconfig_2.0.3       
#>  [69] matrixStats_0.59.0     bitops_1.0-7           evaluate_0.14          lattice_0.20-44       
#>  [73] purrr_0.3.4            bit_4.0.4              tidyselect_1.1.1       magrittr_2.0.1        
#>  [77] R6_2.5.0               IRanges_2.26.0         generics_0.1.0         DBI_1.1.1             
#>  [81] pillar_1.6.1           survival_3.2-11        KEGGREST_1.32.0        scatterplot3d_0.3-41  
#>  [85] RCurl_1.98-1.3         tibble_3.1.2           crayon_1.4.1           utf8_1.2.1            
#>  [89] skmeans_0.2-13         viridis_0.6.1          GetoptLong_1.0.5       blob_1.2.1            
#>  [93] digest_0.6.27          xtable_1.8-4           brew_1.0-6             openssl_1.4.4         
#>  [97] stats4_4.1.0           munsell_0.5.0          viridisLite_0.4.0      askpass_1.1
```




