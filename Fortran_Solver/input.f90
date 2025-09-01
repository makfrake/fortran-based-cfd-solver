subroutine input
use Variabili

implicit none

	ndim=2 !numero dimensioni problema
	neqs=4 ! numero equazioni 

    open(unit=1,file='input.txt')
    read(1,*) 
    read(1,*) 
    read(1,*) mesh_file
    read(1,*) 
    read(1,*) 
    read(1,*) kfinal
    read(1,*) 
    read(1,*) 
    read(1,*) kinf
    read(1,*) 
    read(1,*) 
    read(1,*) kout   
    read(1,*) 
    read(1,*) 
    read(1,*) CFL           
    close(1)
 
 
    open(unit=1,file='inlet.txt')
    read(1,*) 
    read(1,*) 
    read(1,*) ttotin
    read(1,*) 
    read(1,*) 
    read(1,*) ptotin    
    read(1,*) 
    read(1,*) 
    read(1,*) alpha   
    read(1,*) 
    read(1,*) 
    read(1,*) machin
    close(1)
    
    alpha=alpha*(4.*atan(1.0))/(180.)
    
    open(unit=1,file='outlet.txt')
    read(1,*) 
    read(1,*) 
    read(1,*) pexit    
    close(1)    


   



end subroutine




