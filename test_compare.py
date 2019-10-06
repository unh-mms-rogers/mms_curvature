def deepcheck(sidea,sideb):
  try:
    alen = len(a)
  except:
    alen = -1
  try:
    blen = len(b)
  except:
    blen = -1
  if alen != blen:
    return False
  if alen >= 0:
    answer = True
    for i in range(alen):
      try:
        answer &= deepcheck(sidea[i],sideb[i])
      except:
        answer &= sidea == sideb
    return answer
  return sidea == sideb
  
import dataload_parallel as parallel
import dataload_spedas as spedas
import datetime as dt
a=[]
b=[]

timestart = dt.datetime.now()
a.extend(parallel.DataLoad())
timemid = dt.datetime.now()
b.extend(spedas.DataLoad())
timeend = dt.datetime.now()

paralleltime = timemid - timestart
spedastime = timeend - timemid

results = []
matching = True
for i in range(len(a)):
  results.append(deepcheck(a[i],b[i]))
  matching &= results[i]

if matching:
  print("Resulting data matches.")
else:  
  print("Resulting data does not match.  Limited details follow:")
  print(results)

print("parallel runtime: "+str(paralleltime))
print("pyspedas runtime: "+str(spedastime))
