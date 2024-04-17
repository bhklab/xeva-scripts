These files are used for converting XevaSets into CSVs which are used to seed the XevaDB web app database. You will first need to create a xevasets-obj, results, and data folder at the top level of your directory. xevasets-obj will contain the XevaSets to be referenced for conversion, data will contain extra data such as the mutation.xlsx, and the results will be where all of your output will be stored. Note, in the results directory you may need to create subdirectories depending on which script you are running.


R packages needed include Xeva, readxl, reshape2, and data.table
