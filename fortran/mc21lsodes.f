C
C	mc21_lsodes - 2 species, 1 patch type LV master equation solver
C
C   Christopher Klausmeier (klausme1@msu.edu)
C
C	guide to Y's:
C
C	P(N1,N2) -> Y(N1+N2*(N1MAX+1)+1)
C
C
        external F,JAC
        double precision ATOL,RPAR,RTOL,T,TNEXT,TSKIP,TMAX
        integer X,N1,N2,N1MAX,N2MAX,NEQ,MF,NNZ
        double precision, dimension(:), allocatable :: Y,YDOT,RWORK
        integer, dimension(:), allocatable :: IWORK
C
        double precision R1,R2,K1,K2,A12,A21,I1,I2,M1,M2,E ! parameters
C
        common /commonvars/ N1MAX,N2MAX,R1,R2,K1,K2,A12,A21,I1,I2,M1,M2,E
C
        write(*,*) "R1,R2,K1,K2,A12,A21,I1,I2,M1,M2,E"
        read(*,*) R1,R2,K1,K2,A12,A21,I1,I2,M1,M2,E
        write(*,*) "N1MAX,N2MAX,TSKIP,TMAX"
        read(*,*) N1MAX,N2MAX,TSKIP,TMAX
        write(*,*) "RTOL,ATOL,MF,NNZ"
        read(*,*) RTOL,ATOL,MF,NNZ
C
        NEQ=(N1MAX+1)*(N2MAX+1)
        select case (MF)
            case (10)
                LRW=20+16*NEQ
                LIW=30
            case (121)
                LRW=22+3*NNZ+20*NEQ
                LIW=30
            case (222)
                LRW=22+3*NNZ+20*NEQ
                LIW=30
        end select
C
        allocate(Y(NEQ),YDOT(NEQ))
        allocate(RWORK(LRW))
        allocate(IWORK(LIW))
C
        do X=5,10
            RWORK(X)=0.0
            IWORK(X)=0
        enddo
        IWORK(5)=5				! max order
        IWORK(6)=10000000		! maxsteps
C
        ITOL=1                  ! tolerance in scalar form (1)
C
        ITASK=1
        ISTATE=1
        IOPT=1
C
C	input initial conditions
C
        do N2=0,N2MAX
            do N1=0,N1MAX
                read(*,*) Y(N1+N2*(N1MAX+1)+1)
            enddo
        enddo
C
        open(unit=6,file='messages')
        open(unit=7,file='output')
C
        T=0.0
C
C	output data
C
C       write(7,*) "t=",t
        do N2=0,N2MAX
            do N1=0,N1MAX
                write(7,101) Y(N1+N2*(N1MAX+1)+1)
            enddo
        enddo
C
        TNEXT=0.0
        do while (T .LT. TMAX)
            TNEXT=TNEXT+TSKIP
            call DLSODES(F,NEQ,Y,T,TNEXT,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,RWORK,LRW,IWORK,LIW,JAC,MF)
C
C	output data
C
C           write(7,*) "t=",t
            do N2=0,N2MAX
                do N1=0,N1MAX
                    write(7,101) Y(N1+N2*(N1MAX+1)+1)
                enddo
            enddo
            write(6,60) T,IWORK(11),IWORK(12),IWORK(13),IWORK(19),IWORK(20),IWORK(21),IWORK(22)
C
            if (ISTATE .LT. 0) goto 80	!! bad LSODE, no cookie
            T=TNEXT
        enddo
        close(7)
        close(6)
        stop
C
C	error occured
C
 80     write(6,90) ISTATE
        stop
C
C	formats
C
 60     format('t=',F8.1,' steps=',I12,' f-s=',I8,' J-s=',I8,' LU-s=',I8,' Newtons=',I8,'  Newton fails=',I8,' Error test fails=',I8)
 90     format(///' Error halt: ISTATE =',I3)
101     format(F20.10)
C
        end
C
C
C
        subroutine F(NEQ,T,Y,YDOT,RPAR,IPAR)
        double precision RPAR,T,Y,YDOT
        dimension Y(NEQ),YDOT(NEQ)
        integer N1,N2,N1MAX,N2MAX
C
        double precision R1,R2,K1,K2,A12,A21,I1,I2,M1,M2,E
        double precision ITOT1,ITOT2,TOZERO
        double precision INC1,INC2,DEC1,DEC2,EXT
        dimension INC1(0:N1MAX,0:N2MAX),INC2(0:N1MAX,0:N2MAX),DEC1(0:N1MAX,0:N2MAX),DEC2(0:N1MAX,0:N2MAX),EXT(0:N1MAX,0:N2MAX)
C
        common /commonvars/ N1MAX,N2MAX,R1,R2,K1,K2,A12,A21,I1,I2,M1,M2,E
C
C   set up transition rates
C
		TOZERO=0.0
		ITOT1=I1
		ITOT2=I2
        do N2=0,N2MAX
            do N1=0,N1MAX
                DEC1(N1,N2)=(R1*(N1+A12*N2)/K1+M1)*N1*Y(N1+N2*(N1MAX+1)+1)
                DEC2(N1,N2)=(R2*(A21*N1+N2)/K2+M2)*N2*Y(N1+N2*(N1MAX+1)+1)
                EXT(N1,N2)=E*Y(N1+N2*(N1MAX+1)+1)
				TOZERO=TOZERO+EXT(N1,N2)
				ITOT1=ITOT1+M1*N1*Y(N1+N2*(N1MAX+1)+1)
				ITOT2=ITOT2+M2*N2*Y(N1+N2*(N1MAX+1)+1)
            enddo
        enddo
		do N2=0,N2MAX
            do N1=0,N1MAX
                INC1(N1,N2)=(R1*N1+ITOT1)*Y(N1+N2*(N1MAX+1)+1)
                INC2(N1,N2)=(R2*N2+ITOT2)*Y(N1+N2*(N1MAX+1)+1)
            enddo
        enddo
C
C
C   set up YDOTs -- P(N1,N2) = Y(N1+N2*(N1MAX+1)+1)
C
C
C   (N1,N2)=(0,0)
C
        YDOT(1)=-INC1(0,0)-INC2(0,0)+DEC1(1,0)+DEC2(0,1)-EXT(0,0)+TOZERO
C
C   (N1,N2)=(N1,0)
C
        do N1=1,N1MAX-1
            YDOT(N1+1)=INC1(N1-1,0)-INC1(N1,0)-INC2(N1,0)+DEC1(N1+1,0)-DEC1(N1,0)+DEC2(N1,1)-EXT(N1,0)
        enddo
C
C   (N1,N2)=(N1MAX,0)
C
        YDOT(N1MAX+1)=INC1(N1MAX-1,0)-INC2(N1MAX,0)-DEC1(N1MAX,0)+DEC2(N1MAX,1)-EXT(N1MAX,0)
C
C   (N1,N2)=(0,N2)
C
        do N2=1,N2MAX-1
            YDOT(N2*(N1MAX+1)+1)=-INC1(0,N2)+INC2(0,N2-1)-INC2(0,N2)+DEC1(1,N2)+DEC2(0,N2+1)-DEC2(0,N2)-EXT(0,N2)
        enddo
C
C   (N1,N2)=(N1,N2)
C
        do N2=1,N2MAX-1
            do N1=1,N1MAX-1
                YDOT(N1+N2*(N1MAX+1)+1)=INC1(N1-1,N2)-INC1(N1,N2)+INC2(N1,N2-1)-INC2(N1,N2)+DEC1(N1+1,N2)-DEC1(N1,N2)+DEC2(N1,N2+1)-DEC2(N1,N2)-EXT(N1,N2)
            enddo
        enddo
C
C   (N1,N2)=(N1MAX,N2)
C
        do N2=1,N2MAX-1
            YDOT(N1MAX+N2*(N1MAX+1)+1)=INC1(N1MAX-1,N2)+INC2(N1MAX,N2-1)-INC2(N1MAX,N2)-DEC1(N1MAX,N2)+DEC2(N1MAX,N2+1)-DEC2(N1MAX,N2)-EXT(N1MAX,N2)
        enddo
C
C   (N1,N2)=(0,N2MAX)
C
        YDOT(N2MAX*(N1MAX+1)+1)=-INC1(0,N2MAX)+INC2(0,N2MAX-1)+DEC1(1,N2MAX)-DEC2(0,N2MAX)-EXT(0,N2MAX)
C
C   (N1,N2)=(N1,N2MAX)
C
        do N1=1,N1MAX-1
            YDOT(N1+N2MAX*(N1MAX+1)+1)=INC1(N1-1,N2MAX)-INC1(N1,N2MAX)+INC2(N1,N2MAX-1)+DEC1(N1+1,N2MAX)-DEC1(N1,N2MAX)-DEC2(N1,N2MAX)-EXT(N1,N2MAX)
        enddo
C
C   (N1,N2)=(N1MAX,N2MAX)
C
        YDOT(N1MAX+N2MAX*(N1MAX+1)+1)=INC1(N1MAX-1,N2MAX)+INC2(N1MAX,N2MAX-1)-DEC1(N1MAX,N2MAX)-DEC2(N1MAX,N2MAX)-EXT(N1MAX,N2MAX)
C
        return
        end

        subroutine JAC(NEQ,T,Y,J,IA,JA,PDJ)
        double precision T,Y,PDJ
        integer NEQ,J,IA,JA
        integer N1,N2,N1MAX,N2MAX,N1J,N2J
        double precision R1,R2,K1,K2,A12,A21,I1,I2,M1,M2,E
C
        dimension Y(NEQ),PDJ(NEQ),IA(*),JA(*)
C
        common /commonvars/ N1MAX,N2MAX,R1,R2,K1,K2,A12,A21,I1,I2,M1,M2,E
C
C   set up jacobian rates
C
        N1J=mod(J-1,N1MAX+1)
        N2J=(J-1)/(N1MAX+1)
C        write(6,*) J,N1J,N2J
        PDJ(N1J+N2J*(N1MAX+1)+1)=-(R1*N1J+I1+R2*N2J+I2+(R1*(N1J+A12*N2J)/K1+M1)*N1J+(R2*(A21*N1J+N2J)/K2+M2)*N2J+E)
        if (N1J .NE. 0) then
            PDJ((N1J-1)+N2J*(N1MAX+1)+1)=(R1*(N1J+A12*N2J)/K1+M1)*N1J
        end if
        if (N2J .NE. 0) then
            PDJ(N1J+(N2J-1)*(N1MAX+1)+1)=(R2*(A21*N1J+N2J)/K2+M2)*N2J
        end if
        if (N1J .NE. N1MAX) then
            PDJ((N1J+1)+N2J*(N1MAX+1)+1)=R1*N1J+I1
        end if
        if (N2J .NE. N2MAX) then
            PDJ(N1J+(N2J+1)*(N1MAX+1)+1)=R2*N2J+I2
        end if
        PDJ(1)=PDJ(1)+E
        return
        end
