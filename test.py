__author__ = 'Tiannan Guo, ETH Zurich 2015'


a = {}
a[1]=2
a[2]=2
a[3]=4
a[4]=5
a[5]=6
a[6]=[1,2,3,4,5]
a[7]=[5,6,7]
b=a
c=a.copy()

del b[1]
