Program temp

  Implicit None
  !Variables
  Integer:: i,n,Nfin
  Integer, parameter::Np=100000
  Real*8,Dimension(Np)::Tab_Ep,sigma
  Real*8,Dimension(100)::Hist
  Real*8,Dimension(:),Allocatable::Tab_Tint
  Real*8::tau,dt,T,R,u,v,a,b,c,cst
  Real*8,Parameter::pi=3.1415927
  real :: start, finish

  tau=5.
  dt=1./100
  T=300.
  R=280.
  Nfin=floor(5*tau/dt)
  Allocate(Tab_Tint(Nfin+1))

  Open(1,file="Tint.dat",form="formatted",status="unknown")

  !debut timer 
  call cpu_time(start)


	     Tab_Ep=0.
  !tirages au sort : 500000 Ep entre 100000J et 200000J
     Do i=1,Np
        Call Random_number(Tab_Ep(i))
        Tab_Ep(i)=Tab_Ep(i)*100000.+(1-Tab_Ep(i))*200000.
     End Do
     !Calcul Tint
     Tab_Tint=0.
     Do i=1,Np
        Tab_Tint(1)=Tab_Ep(i)+Tab_Tint(1)
     End Do
     Tab_Tint(1)=Tab_Tint(1)/Np
     Write(1,'(1I5,a,1F9.2)')1,' ',Tab_Tint(1)
     !constantes
     a=2*pi
     b=1./(1+2*dt/tau)
     c=dt/tau*R*T
     
     !calcul des sigmas
 

     !Entre t et t+dt
     Do n=2,Nfin+1
          Do i=1,Np/2
           Call Random_number(u)
           Call Random_number(v)
           cst=SQRT(-2.*Log(1-u))
           sigma(i)=cst*Cos(a*v)
           sigma(i+Np/2)=cst*sin(a*v)
        enddo
        Do i=1,Np
           Tab_Ep(i)=b*(Tab_Ep(i)+c*(1+sigma(i)**2)+2*sigma(i)*SQRT(c*Tab_Ep(i)))
           Tab_Tint(n)=Tab_Tint(n)+Tab_Ep(i)
        End Do
        !Tab_Tint(n)=Tab_Tint(n)/Np
     End Do
     Do n=2,Nfin+1
        Write(1,'(1I5,a,1F9.2)')n,' ',Tab_Tint(n)/Np
     end Do
     Close(1)
  !end timer
  call cpu_time(finish)
   print '("Time = ",f6.3," seconds.")',finish-start
  !Création d'un histogramme

  Hist=0.
  Open(2,file="hist.dat",status="unknown")
  Do i=1,Np
     if (int(Tab_Ep(i))<5*R*T) then
        Hist(int(int(Tab_Ep(i))*100/(5*R*T)+1))=Hist(int(int(Tab_Ep(i))*100/(5*R*T)+1))+1
        !print*,int(int(Tab_Ep(i))*100/(5*R*T)+1)
     end if
  end do
  Do i=1,100
     Write(2,*)i,Hist(i)
  end do
  Close(2)

  Deallocate(Tab_Tint)

End Program temp
