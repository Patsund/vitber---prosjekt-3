import oppgave1 as o1
import oppgave2 as o2
import oppgave3 as o3

def printOptions():
	outStr = "\n1) Oppgave 1"
	outStr+= "\n2) Oppgave 2"
	outStr+= "\n3) Oppgave 3"
	outStr+= "\n0) Avslutt"
	outStr+= "\nDitt valg:"
	return outStr 

if __name__ == "__main__":
	switch = -1
	while(switch):
		try:
			switch = int(input(printOptions()))
		except:
			pass
		if switch == 1:
			o1.oppgave1()

		elif switch == 2:
			o2.oppgave2()
		elif switch == 3:
			o3.oppgave3()
		elif switch == 0:
			break
		else:
			print("Invalid input\n")