import time

fid = open('/home/hyde/testcron.txt', 'w');

fid.write('time is '+time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())+'\n')

