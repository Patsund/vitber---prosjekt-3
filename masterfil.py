import oppgave1 as o1
import oppgave2 as o2
import Oppgave3 as o3

def printOptions():
	outStr = "\n1) Oppgave 1"
	outStr+= "\n2) Oppgave 2"
	outStr+= "\n3) Oppgave 3"
	outStr+= "\n0) Avslutt"
	outStr+= "\nDitt valg:"
	return outStr

if __name__ == "__main__":
	switch = -1
	save = -1
	while(switch):
		try:
			switch = int(input(printOptions()))
			save = int(input("\n1) Lagre figurer \n2) Vis figurer \n 0) Avslutt"))
		except:
			pass
		if switch == 1:
			if save == 2:
				o1.oppgave1()
			elif save == 3:
				o1.oppgave1(True)
			else:
				print("Du oppga ikke en gyldig verdi for å lagre/vise figurer")
				pass
		elif switch == 2:
			o2.oppgave2()
		elif switch == 3:
			if save==2:
				o3.oppgave3()
			elif save==1:
				o3.oppgave3(True)
			else:
				print("Du oppga ikke en gyldig verdi for å lagre/vise figurer")
				pass
		elif switch == 0 or save == 0:
			break
		else:
			print("Invalid input\n")