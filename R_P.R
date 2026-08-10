#1.安装reticulate包和moniconda
install.packages("reticulate")
library("reticulate")
#install_miniconda()

#2.加载conda虚拟环境或python
#使用conda虚拟环境
use_condaenv("r-reticulate")
#使用已安装的python
#use_python("C:\\Users\\38961\\AppData\\Local\\Programs\\Python\\Python310\\python.exe")
##也可以通过Options点击设置

#3.在R环境下安装python包!!!!
py_install("scrublet")
py_install("pip")
py_install("fa2")
#查看已经安装的包
py_list_packages()
#4.在R中执行python文件
soure_python("found.py")
