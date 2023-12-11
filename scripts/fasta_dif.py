#list_1=[]
#file = open('fff.txt')
#file2 = open('ff2.txt')
#for line in file:
#	 count=count+1	


#python fasta_dif.py Or_aa.fasta Or_nn_v2.fasta    
import sys
input2=sys.argv[1]
input1=sys.argv[2]
count=0
err=0
name_aa=""
with open(input1,'r+') as f1:
	for line1 in f1:
		
		count=count+1
		if count%2==1:
			name_nn=line1
			with  open(input2) as f4 :
				for  line4 in  f4:
					if name_nn==name_aa:
						len2=len(line4)-1
						break
					name_aa=line4
		else :
			len1=len(line1)-1
			if len2*3!=len1:
				err=1
				print("length different")
				print(name_nn + " length: " + str(len2) + "*3=" + str(len2*3))
				print(name_aa + " length: " + str(len1))
				print("-----------------------------")
				print("\n")

	#	if count%2==0:
	#		len1=len(line1)-1
	#		len2=len(line2)-1
	#		if len2*3!=len1:
	#			print("nooooooooooooo")
	#			print line1
	#			print line2
	#			print len1
	#			print len2*3
if err==1:
	sys.exit()