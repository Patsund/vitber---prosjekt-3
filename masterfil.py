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

def printSave():
	outStr = "\n1) Vis figurer"
	outStr+= "\n2) Lagre figurer"
	outStr+= "\n0) Tilbake"
	outStr+= "\nDitt valg:"
	return outStr

if __name__ == "__main__":
	switch = -1
	save = -1
	while(switch):
		try:
			switch = int(input(printOptions()))
			save = int(input(printSave()))
		except:
			pass
		if save==0:
			pass
		if switch == 1:
			if save == 1:
				o1.oppgave1()
			elif save == 2:
				o1.oppgave1(True)
			else:
				print("Du oppga ikke en gyldig verdi for å lagre/vise figurer")
				pass
		elif switch == 2:
			datapath = 'C:/Users/Patrik/Downloads/NorKyst-800m.nc'
			if save == 1:
				o2.oppgave2(datapath)
			elif save == 2:
				o2.oppgave2(datapath, True)
			else:
				print("Du oppga ikke en gyldig verdi for å lagre/vise figurer")
				pass
		elif switch == 3:
			datapath = 'C:/Users/Patrik/Downloads/NorKyst-800m.nc'
			if save == 1:
				o3.oppgave3(datapath)
			elif save == 2:
				o3.oppgave3(datapath, True)
			else:
				print("Du oppga ikke en gyldig verdi for å lagre/vise figurer")
				pass
		elif switch == 0:
			break
		else:
			print("Du oppga ikke en gyldig verdi for å velge oppgave\n")