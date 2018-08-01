#! /usr/bin/python3
import matplotlib.pyplot as plt

f=open('rtfile.txt')
#f2=open('x.txt','w')
#f3=open('z.txt','w')

x,y,z = [], [], []
for l in f:
    row = l.split()
    x.append(float(row[0]))
    z.append(float(row[1]))
    y.append(float(row[2]))

#for item in x:
# f2.write("%s\n" % item)

#for item in y:
# f3.write("%s\n" % item)



plt.plot(x,y)

plt.show()