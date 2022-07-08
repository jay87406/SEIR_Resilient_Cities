# 环境

pip install flask

pip install networkx

pip install numpy

pip install scipy

# 内容

透过server.py来启动，大概会运行1-3 mins

其他的文件都是用来调用的(有些不重要，但那时候想说急着看所以全部打包了)

## 输入

目前需要输入的参数有5个，其他都用内置参数

| 变量名                   | 类型  | 示例 | 描述                   |
| ------------------------ | ----- | ---- | ---------------------- |
| MEAN_INTRACOHORT_DEGREE  | int | 18   | 社群网络内节点degree   |
| PCT_CONTACTS_INTERCOHORT | float | 0.20 | 社群网络接触后感染可能 |
| INIT_EXPOSED             | int | 30   | 初始感染节点数         |
| R0                       | float | 2.5  | 描述感染程度的值       |
| T                        | int   | 30   | 运行时间(day)          |

里面用之前的方式写好了

```python
	MEAN_INTRACOHORT_DEGREE = float(request_body['MEAN_INTRACOHORT_DEGREE']) or 18
    PCT_CONTACTS_INTERCOHORT = float(request_body['PCT_CONTACTS_INTERCOHORT']) or 0.20
    INIT_EXPOSED = float(request_body['INIT_EXPOSED']) or 30
    #------------------------------------------
    R0 = float(request_body['R0']) or 2.5
    T = int(request_body['T']) or 30
```

