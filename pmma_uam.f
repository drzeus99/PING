	program pmma_uam
	implicit none
C Ability to set up PMMA chains around a hollow sphere in a simulation box of size lx,ly,lz
C lx,ly,lz depend on the %alumina mixed into the pmma. The size of the hollow sphere will 
C always be 1 angstrom diameter. But in this version the hollow sphere will be skipped.
	integer chainnum, atomnum, unitnum, kuhncount
	parameter (chainnum=30, unitnum=100, atomnum=7,kuhncount=7)
c Chainnum is the number of chains, unitnum is the number of pmma units in the chain
c atomnumber is the number of particles in a unit	
	integer bondnum, anglenum, dihedralnum, impropernum
	integer atomtypenum, bondtypenum, angletypenum, dihedraltypenum, impropertypenum
	parameter (atomtypenum=6, bondtypenum=6, angletypenum=6, dihedraltypenum=6,impropertypenum=1)
	double precision atoms(5,atomnum), adition(3), masses(atomtypenum)
	double precision chainatoms(3,unitnum*atomnum), pbcchain(3,unitnum*atomnum)
	integer imageflag(3,unitnum*atomnum)
	double precision bondcoef(2,bondtypenum), anglecoef(2,angletypenum)
	double precision dihedralcoef(5,dihedraltypenum),impropercoef(2,impropertypenum)
	integer bondinit(3,6), bonds(3,7), angleinit(4,7), angles(4,11)
	integer dihedralinit(5,6), dihedrals(5,13), impropers(5)
	double precision hspherer, latticesize
	parameter (hspherer=15.0, latticesize=75.845720286)
C hsperer, latticesize and minshift are in angstrom.
C hspher is the radius of the hollow sphere in the center of the simulation box	
	double precision shiftdist, shiftunit(3), shift(3), dist2
	integer i, j, k, l,n, sphereflag, count, unitcount
	data atoms/3, .11, 0, 0, 0,
     &	2, 0, -1.52429914864554, 0,	.2193447182826,
     &	1, 0, .333005967244697,	-1.23820341059349,	-.8529708903439,
     &	4, .31, .75,	0,	1.29903810567666,
     &	5, -.37, .262586755537238,	.528275496308508,	2.28481360198517,
     &	6, -.31, 2.03631224555881,	-.619208163705874,	1.38199626938253,
     &	1, .26, 2.43869486125905,	-.42813272480992,	2.74084814035043/,
     &	adition/.333005967244697,	1.23820341059349,	-.852970890343899/,
     &	masses/12.0107, 12.0107, 12.0107, 12.0107, 15.9994, 15.9994/,
     &	bondinit/1, 1, 3,
     &	2, 1, 2,
     &	3, 1, 4,
     &	4, 4, 5,
     &	5, 4, 6,
     &	6, 6, 7/,
     &	bonds/2, 1, -5,
     &	1, 1, 3,
     &	2, 1, 2,
     &	3, 1, 4,
     &	4, 4, 5,
     &	5, 4, 6,
     &	6, 6, 7/,
     &	bondcoef/368, 1.539,
     &	300, 1.549, 
     &	326, 1.517,
     &	968, 1.209,
     &	471, 1.360,
     &	342, 1.446/,
     &	angleinit/2, 2, 1, 3,
     &	2, 3, 1, 4,
     &	2, 2, 1, 4,
     &	4, 1, 4, 5,
     &	3, 1, 4, 6,
     &	5, 6, 4, 5,
     &	6, 7, 6, 4/,
     &	angles/1, -6, -5, 1,
     &	2, -5, 1, 4,
     &	2, -5, 1, 3,
     &	2, 2, 1, 3,
     &	2, 2, 1, 4,
     &	2, 2, 1, -5,
     &	2, 3, 1, 4,
     &	3, 6, 4, 1,
     &	4, 5, 4, 1,
     &	5, 6, 4, 5,
     &	6, 7, 6, 4/,
     &	anglecoef/89.5, 113.3,
     &	87.9, 109.47,
     &	74.5, 111.4,
     &	63.3, 125.6,
     &	126.5, 123.0,
     &	84.8, 116.4/,
     &	dihedralinit/4, 3, 1, 4, 5,
     &	4, 2, 1, 4, 5,
     &	3, 3, 1, 4, 6, 
     &	3, 2, 1, 4, 6,
     &	5, 1, 4, 6, 7,
     &	6, 7, 6, 4, 5/,
     &	dihedrals/1, 1, -5, -6, -4,
     &	2, 1, -5, -6, -3,
     &	1, 2, 1, -5, -6,
     &	1, 3, 1, -5, -6,
     &	2, 4, 1, -5, -6,
     &	4, 5, 4, 1, -5,
     &	3, 6, 4, 1, -5,
     &	4, 3, 1, 4, 5,
     &	4, 2, 1, 4, 5,
     &	3, 3, 1, 4, 6,
     &	3, 2, 1, 4, 6,
     &	6, 7, 6, 4, 5,
     &	5, 7, 6, 4, 1/,
     &	dihedralcoef/.00043047, .88728, -.0058908, 3.48522, .0083951,
     &	-.69938, .88749, 1.39251, 3.48489, .010033,
     &	-.00011517, 1.08001, .00042681, 1.76002, -.00036242,
     &	1.90050, -.78234, -3.80681, 1.04477, .0094108,
     &	2.22989, 4.54992, -4.45939, -.63999, -.00051598,
     &	2.20044, 1.04739, -4.40604, -1.39492, .0085622/,
     &	impropers/1, 4, 1, 5, 6/,
     &	impropercoef/0.735006480784105,	0/

	
C write the header portion of the data file
	open (25,file='data.pmma90')
	write (25,*) 'this file contains the pmma polymer chains in the UAM format'
C Sets up the simulation box size	
	write (25,*) -latticesize/2d0, latticesize/2d0, ' xlo xhi' 
	write (25,*) -latticesize/2d0, latticesize/2d0, ' ylo yhi'
	write (25,*) -latticesize/2d0, latticesize/2d0, ' zlo zhi'
C Sets up the number of atoms
	write (25,*) chainnum*unitnum*atomnum, ' atoms'
C sets up the number of bonds, angles, dihedrals and impropers
	bondnum=chainnum*(unitnum-1)*atomnum+chainnum*(atomnum-1)
	anglenum=chainnum*(unitnum-1)*(atomnum+4)+chainnum*(atomnum)
	dihedralnum=chainnum*(atomnum-1)+chainnum*(unitnum-1)*(2*atomnum-1)
	impropernum=chainnum*unitnum
	write (25,*) bondnum, ' bonds'
	write (25,*) anglenum, ' angles'
	write (25,*) dihedralnum, ' dihedrals'
	write (25,*) impropernum, ' impropers'
C set up the the number of differnt types of bonds, angles, dihedrals and impropers
	write (25,*) atomtypenum, ' atom types'
	write (25,*) bondtypenum, ' bond types'
	write (25,*) angletypenum, ' angle types'
	write (25,*) dihedraltypenum, ' dihedral types'
	write (25,*) impropertypenum, ' improper types'
	
	do 2 i=1,3
		write (25,*)
 2	continue
	
	call srand(14)
c	call srand(7)
c	call srand(3)
	write (25,*) 'Atoms'
	write (25,*)
	open (30,file='pmma90.xyz')
C write the xyz file header
	write (30,*) chainnum*unitnum*atomnum 
	write (30,*) 'atoms'
C atom counter for polymer data file
	count=1
C build the atoms in the box
	do 15 i=1,chainnum
c place in first unit of the new chain randomly 	
 1		do 16 j=1,3
			chainatoms(j,1)=latticesize*rand(0)-latticesize/2.0
			call pbc(chainatoms(j,1),latticesize,pbcchain(j,1), imageflag(j,1))
 16		continue
		if (hspherer .ne. 0) then
			dist2=pbcchain(1,1)**2+pbcchain(2,1)**2+pbcchain(3,1)**2
			if (dist2 .le. hspherer**2) goto 1
		endif
				
 		do 20 j=2,atomnum
C just do the adition of atoms to chainatoms and pbcchain	
			do 22 k=1,3
				chainatoms(k,j)=chainatoms(k,1)+atoms(k+2,j)
				call pbc(chainatoms(k,j),latticesize,pbcchain(k,j),imageflag(k,j))
 22			continue
C checking to make sure atoms are not in the hollow sphere 
			if (hspherer .ne. 0) then
				dist2=pbcchain(1,j)**2+pbcchain(2,j)**2+pbcchain(3,j)**2
				if (dist2 .le. hspherer**2) goto 1
			endif	
 20		continue	
		unitcount=1
			
		do 10 j=2,unitnum
			sphereflag=0
			do 5 k=1,atomnum
c just do the adition of atoms to chainatoms and pbcchain
				dist2=0
				do 7 l=1,3
					if (k==1) then
						chainatoms(l,(j-1)*atomnum+k)=chainatoms(l,(j-1)*atomnum-5)+adition(l)
					else 
						chainatoms(l,(j-1)*atomnum+k)=chainatoms(l,(j-1)*atomnum+1)+atoms(l+2,k)						
					endif
					call pbc(chainatoms(l,(j-1)*atomnum+k),latticesize,pbcchain(l,(j-1)*atomnum+k), imageflag(l,(j-1)*atomnum+k))
					if  (sphereflag .eq. 1) goto 7				
C above if statement skips code below					
					if (hspherer .ne. 0) then
						dist2=pbcchain(l,(j-1)*atomnum+k)**2+dist2
					endif
 7				continue
				if (sphereflag .eq. 1) goto 5				
C above if statement skips code below					
				if (hspherer .ne. 0) then
					if (dist2 .le. hspherer**2) sphereflag=1
				endif									
c check to see if atoms are in the sphere
c if so trigger sphereflag												
 5			continue
 
			if (sphereflag .eq. 1) then
				call mirror(unitnum,atomnum,j-1,chainatoms, pbcchain, atoms, adition, latticesize, imageflag)
				unitcount=0
			elseif (unitcount .eq. kuhncount) then
				call mirror(unitnum,atomnum,j-1,chainatoms, pbcchain, atoms, adition, latticesize, imageflag)
				unitcount=0
			else
				unitcount=unitcount+1
			endif
 10		continue
 
		if (hspherer .ne. 0) then
c checks to make sure no atoms of the current chain are in the sphere
			do 18 j=2,unitnum
c Because of how the chain is built the 1st unit cannot be in the sphere so start with unit 2
				do 78 k=1,atomnum
c Check the atoms in each unit
					dist2=0
					do 12 l=1,3
						dist2=dist2+pbcchain(l,(j-1)*atomnum+k)**2
 12					continue
					if (dist2 .le. hspherer**2) then
C calculate the shift unit vector						
						dist2=sqrt(dist2)
						do 9 l=1,3
							shiftunit(l)=pbcchain(l,(j-1)*atomnum+k)/dist2
 9						continue
						shiftdist=hspherer+.1
C calculate the chainatoms, imageflag and pbcchain
						do 13 l=1,3
							chainatoms(l,(j-1)*atomnum+k)=chainatoms(l,(j-1)*atomnum+k)+shiftdist*shiftunit(l)-pbcchain(l,(j-1)*atomnum+k)
							call pbc(chainatoms(l,(j-1)*atomnum+k), latticesize, pbcchain(l,(j-1)*atomnum+k), imageflag(l,(j-1)*atomnum+k))
 13						continue
					endif
 78				continue
 18			continue
		endif
				
C write the final chain structure in the data file and in a .xyz file
		do 6 j=1,unitnum		
			do 3 k=1,atomnum
C data file			
				write (25,*) count, i, int(atoms(1,k)), atoms(2,k), pbcchain(1,(j-1)*atomnum+k),
     & 				pbcchain(2,(j-1)*atomnum+k), pbcchain(3,(j-1)*atomnum+k), imageflag(1,(j-1)*atomnum+k),
     &				imageflag(2,(j-1)*atomnum+k), imageflag(3,(j-1)*atomnum+k)
				count=count+1
C xyz file	
				if	(atoms(1,k) .lt. 5) then
					write (30,*) 6, pbcchain(1,(j-1)*atomnum+k), pbcchain(2,(j-1)*atomnum+k),
     &					pbcchain(3,(j-1)*atomnum+k)
				else
					write (30,*) 8, pbcchain(1,(j-1)*atomnum+k), pbcchain(2,(j-1)*atomnum+k),
     & 					pbcchain(3,(j-1)*atomnum+k)
				endif
 3			continue
 6		continue
C Reset atoms and adition
C        call reset(atoms, adition, atomnum)
 15	continue
	close (30)
	write (25,*)
	
c Writing the masses
	write (25,*) 'Masses'
	write (25,*)
	do 25 i=1,atomtypenum
		write (25,*) i, masses(i)
 25	continue
	write (25,*)
 
C write the bonds  secton
	write (25,*) 'Bonds'
	write (25,*)
	count=1
	do 35 i=1,chainnum
		do 36 j=1,atomnum-1
			write (25,*) count, bondinit(1,j), bondinit(2,j)+(i-1)*unitnum*atomnum, 
     &			bondinit(3,j)+(i-1)*unitnum*atomnum
			count=count+1
 36		continue
		do 37 j=2,unitnum
			do 39 k=1,atomnum
				write (25,*) count, bonds(1,k), bonds(2,k)+(i-1)*unitnum*atomnum+(j-1)*atomnum,
     & 				bonds(3,k)+(i-1)*unitnum*atomnum+(j-1)*atomnum
				count=count+1
 39			continue
 37		continue
 35	continue
	write (25,*)

C write the bond coeficients section
	write (25,*) 'Bond Coeffs'
	write (25,*)
	do 30 i=1,bondtypenum
		write (25,*) i, bondcoef(1,i), bondcoef(2,i)
 30	continue
	write (25,*)

C write the angles section 
	write (25,*) 'Angles'
	write (25,*)
	count=1
	do 40 i=1,chainnum
		do 41 j=1,atomnum
			write (25,*) count, angleinit(1,j), angleinit(2,j)+(i-1)*unitnum*atomnum,
     & 			angleinit(3,j)+(i-1)*unitnum*atomnum, angleinit(4,j)+(i-1)*unitnum*atomnum
			count=count+1
 41		continue
		do 42 j=2,unitnum
			do 44 k=1,atomnum+4
				write (25,*) count, angles(1,k), angles(2,k)+(i-1)*unitnum*atomnum+(j-1)*atomnum,
     &				angles(3,k)+(i-1)*unitnum*atomnum+(j-1)*atomnum,
     & 				angles(4,k)+(i-1)*unitnum*atomnum+(j-1)*atomnum
				count=count+1
 44			continue
 42		continue
 40	continue
	write (25,*)

C write the angle coeficients section 
	write (25,*) 'Angle Coeffs'
	write (25,*)
	do 43 i=1,angletypenum
		write (25,*) i, anglecoef(1,i), anglecoef(2,i)
 43	continue
	write (25,*)
	
C write the dihedral section
	write (25,*) 'Dihedrals'
	write (25,*)
	count=1
	do 45 i=1,chainnum
		do 46 j=1,atomnum-1
			write (25,*) count, dihedralinit(1,j), dihedralinit(2,j)+(i-1)*unitnum*atomnum,
     & 			dihedralinit(3,j)+(i-1)*unitnum*atomnum, dihedralinit(4,j)+(i-1)*unitnum*atomnum,
     & 			dihedralinit(5,j)+(i-1)*unitnum*atomnum
			count=count+1
 46		continue
		do 47 j=2,unitnum
			do 48 k=1,2*atomnum-1
				write (25,*) count, dihedrals(1,k),
     & 				dihedrals(2,k)+(i-1)*unitnum*atomnum+(j-1)*atomnum, 
     &				dihedrals(3,k)+(i-1)*unitnum*atomnum+(j-1)*atomnum, 
     &				dihedrals(4,k)+(i-1)*unitnum*atomnum+(j-1)*atomnum,
     & 				dihedrals(5,k)+(i-1)*unitnum*atomnum+(j-1)*atomnum
				count=count+1
 48			continue
 47		continue
 45	continue
	write (25,*)


C write the dihedral coefficient section
	write (25,*) 'Dihedral Coeffs'
	write (25,*)
	do 49 i=1,dihedraltypenum
		write (25,*) i, dihedralcoef(1,i), dihedralcoef(2,i), dihedralcoef(3,i), dihedralcoef(4,i),
     & 		dihedralcoef(5,i)
 49	continue
	write (25,*)
	
	
C write impropers section
	write (25,*) 'Impropers'	
	write (25,*)
	count=1
	do 53 i=1,chainnum
		do 52 j=1,unitnum
			write (25,*) count, impropers(1), 
     &			impropers(2)+(i-1)*unitnum*atomnum+(j-1)*atomnum,
     &			impropers(3)+(i-1)*unitnum*atomnum+(j-1)*atomnum,
     &			impropers(4)+(i-1)*unitnum*atomnum+(j-1)*atomnum,
     &			impropers(5)+(i-1)*unitnum*atomnum+(j-1)*atomnum
 			count=count+1
 52		continue
 53	continue
	write (25,*)	

C write the improper coefficient section
	write (25,*) 'Improper Coeffs'
	write (25,*)
	do 51 i=1,impropertypenum
		write (25,*) i, impropercoef(1,i), impropercoef(2,i)
 51	continue
	write (25,*)

	close (25)
	end

	subroutine mirror(unitnum,atomnum,j,chainatoms, pbcchain, atoms, adition,latticesize, imageflag)
	implicit none
C global variables
	integer unitnum, atomnum, j
C j is the current unitnumber	
	double precision chainatoms(3,unitnum*atomnum), pbcchain(3,unitnum*atomnum)
	double precision atoms(5,atomnum), adition(3), latticesize
	integer	imageflag(3,unitnum*atomnum)
C local variables
	double precision oldunit(3,atomnum), newunit(3,atomnum)
C oldunit and newunit are in terms of the pbc positions
C old is the current unit in the pbc and new is the transformed one.
	double precision mirrorpoint(3), dist2, dot
	double precision atomtemp(3), tempatom(3)
	integer i,k
C this algorithm calculates a reflection when the plane of the mirror is parralel to a sphere defined
C by the vector mirrorpoint. mirrorpoint is also the normal vector of the mirror. 

C Calculate dist2
	dist2=0
	do 55 i=1,3
		mirrorpoint(i)=chainatoms(i,(j*atomnum)-5)
C		write (*,*) 'the mirrorpoint at index ', i, ' is ', mirrorpoint(i)
		dist2=dist2+mirrorpoint(i)**2
 55	continue

C Extract oldunit from pbcchain
C atomtemp is in oldunit coordinate system
C tempatom is in newunit coordinate system
C calculate newunit using mirror forumla
C tempatom=atomtemp - 2 * (atomtemp (dot) mirrorpoint - dist2) / dist2 * mirrorpoint
c Use the newunit and oldunit to calculate the new position of chainatoms
C Than calculate the pbcchain using the new chainatoms
	do 59 i=1,atomnum
		do 60 k=1,3
			oldunit(k,i)=chainatoms(k,i+j*atomnum)
			atomtemp(k)=oldunit(k,i)
C			write (*,*) 'the atom at position ', i, ' index ', k, ' is ', atomtemp(k)
 60		continue
C Using mirror formula
        dot=0.0
		do 68 k=1,3
			dot=dot+atomtemp(k)*mirrorpoint(k)
68		continue
		do 65 k=1,3
            tempatom(k)=atomtemp(k)-2*(dot-dist2)/dist2*mirrorpoint(k)
C            write (*,*) 'the mirrored atom at position ', i, ' index '
C     &			, k, ' is ',tempatom(k)
 65		continue
		do 70 k=1,3
			newunit(k,i)=tempatom(k)
c			chainatoms(k,i+j*atomnum)=chainatoms(k,i+j*atomnum)+newunit(k,i)-oldunit(k,i)
			chainatoms(k,i+j*atomnum)=newunit(k,i)
			call pbc(chainatoms(k,i+j*atomnum),latticesize,pbcchain(k,i+j*atomnum), imageflag(k,i+j*atomnum))
 70		continue
 59	continue
 
C Calculate adition
	do 80 k=1,3
		adition(k)=newunit(k,1)-mirrorpoint(k)
C remember mirrorpoint is where the pmma unit is attached to the old chain
 80	continue
	
C Calculate atoms (The 1st column x,y,z values are all 0,0,0)
	do 90 i=2,atomnum
		do 100 k=1,3
			atoms(k+2,i)=newunit(k,i)-newunit(k,1)				
 100	continue
 90	continue
	return
	end

	subroutine pbc(atompos,latticesize, pbcpos, imageflag)
	implicit none
C global variables
	double precision atompos, latticesize, pbcpos
	integer imageflag
C local variables
	double precision min, max, delta
C add image flag into this code
	max=latticesize/2.0
	min=-latticesize/2.0
	delta=max-min
	
	if (atompos .lt. min) then	
		pbcpos=atompos+delta+delta*int(abs((atompos-min)/delta))
        imageflag=-1-int(abs((atompos-min)/delta))
	elseif (atompos .ge. max) then
		pbcpos=atompos-delta-delta*int(abs((atompos-max)/delta))
        imageflag=1+int(abs((atompos-max)/delta))
	else
		pbcpos=atompos
        imageflag=0
	endif
	return
	end

	subroutine reset(atoms, adition, atomnum)
C global variables
	integer atomnum
	double precision atoms(5,atomnum), adition(3)
C only need to reset the position data of atoms so only need to reset rows 3-5

	atoms(3,1)=0
	atoms(4,1)=0
	atoms(5,1)=0
	
	atoms(3,2)=-1.52429914864554
	atoms(4,2)=0
	atoms(5,2)=.2193447182826
	
	atoms(3,3)=.333005967244697
	atoms(4,3)=-1.23820341059349
	atoms(5,3)=-.8529708903439

	atoms(3,4)=.75
	atoms(4,4)=0
	atoms(5,4)=1.29903810567666
	
	atoms(3,5)=.262586755537238
	atoms(4,5)=.528275496308508
	atoms(5,5)=2.28481360198517
	
	atoms(3,6)=2.03631224555881
	atoms(4,6)=-.619208163705874
	atoms(5,6)=1.38199626938253
	
	atoms(3,7)=2.43869486125905
	atoms(4,7)=-.42813272480992
	atoms(5,7)=2.74084814035043

	adition(1)=.333005967244697
	adition(2)=1.23820341059349
	adition(3)=-.852970890343899
	return
	end
