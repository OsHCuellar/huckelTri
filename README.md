

# huckelTri


        Osvaldo Hernandez Cuellar


FORTRAN program that computes and diagonalize the Hückel Hamiltonian for some molecules with with π delocalized electrons and python codes for its analysis. 

## Code structure

```bash
/huckelTri/
    3triangulene.inp  
    Global.inp
    huckel
    inputExamples
        0explicit.inp  
        1nanotube.inp  
        2polyene.inp  
        3triangulene.inp  
        /analysis/
            0coefflvl.py  
            visualization3.py
        Global.inp
    makefile
    output/
        Eigenvectors/
        0explicit.inp
        EnergyLvl.dat
        Global.inp
    README
    src/
        dsyev.f90  
        jacobiMethod.f90  
        main.f90  
        math.f90  
        objs/  
            dsyevdiag.mod  
            jacobidiag.mod  
            main.o    
            math.o         
            start.mod  
            types.o       
            writting.o
            dsyev.o        
            jacobiMethod.o  
            math.mod  
            readMakeInp.o  
            types.mod  
            writting.mod
        readMakeInp.f90  
        types.f90  
        writting.f90    
```

## Usage

### Compile

Because the program uses a subrotuine from the lapack, you have to modify the FFLAGS2 on the makefile to link correclty
the lapack library on your computer with this program.

To install lapack and libblas:
```bash
sudo apt-get install libblas-dev liblapack-dev
```

After doing this you just have to exectute the makefile by writting:
```bash
make 
```
on the directory `/huckelTri/` . This will create the executable `huckelTri` on the same directory.

### Input

The program can create and diagonalize the Hückel Hamiltonian for: 

    A) Explicity defined system
    B) one of the three predefine structures:
        1) nanotube 
        2) polyene
        3) triangulene 

To choose which of the system it is going to solve you have to modify the file `Global.inp`, where you can choose the
type of system, if you want to debug the program, and the diagonalization method.

After chossing the method, you have to verify and modify the input file for the system you want to work with. The input 
file for each system listed above ()in the same order are:

    A) 0explicit.inp
    B)
        1) 1nanotube.inp
        2) 2polyene.inp
        3) 3triangulene.inp

                *** NOTE: all the input files contains the explanations and instructions on 
                    how to modify them, and what does each variable represents.

### Output

In the `/huckelTri/output/` directory are going to be all the outpit files as well as the input files used to generate
them. The EnergyLvl.dat file contains all the energy levels in ascending order, while inside the
`huckelTri/output/Eigenvectors/` directory, there is one file for each energy level called
`/huckelTri/output/Rigenvectors/EigenvectorsN.dat`, where the N indicates the energy level to which it corresponds. Inside
the `EigenvectorsN.dat`, there are all the coefficients of the spin wavefunction for the corresponding energy level.

Regarding the input files created in the `/huckel/output/` directory, besides the `Global.inp` and the corresponding
input file for the system, there is always going to be a `0explicit.inp` file. This file is created automatically and it
correspond to the system specified before, whether you specify it explictly or not. This file can be helpful sometimes 
because it contains information about your system that might not be specified before, like for example the total number 
of atoms, the total number of bonds, among others.

### Execute

To execute the program you just have to execute `./huckelTri` program created after compiling, inside the `/huckelTri/`
directory. To change the system, you just have to modify the input files, recompiling is not necesarry. 

### Analize results

As an extra, to analize some of the results for the triangulene, inside the `/huckelTri/analysis/` directory, there are
some shell scripts and python program to help the analysis of the result. 
    0coefflvl.py
        python code that creates the `0coefflvl.dat` file which contains the a list of the energy levels with the total
        number of coefficients equal to zero for the spin wavefunction for the given energy level. To exectute it it has
        to be in the directory of the Eigenvectors that contains all the eigenvctors and int eh previous directory it
        has to be the corresponfing `oexplicit.dat` file.
    visualization3.py
        python code that makes a visual representation of the spin wave function of the trinalgule specified inside it, 
        on the variable "eigenVec".

            *** NOTE: To execute this programs yopu have ot verify that you have installed 
                the needed pyhton libraries.



