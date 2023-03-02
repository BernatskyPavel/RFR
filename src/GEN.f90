module GEN
   USE RFR, ONLY: GENER, BUCK

   PUBLIC :: VARGEN
contains
! АЛГОРИТМ ГЕНЕРИРОВАНИЯ СКАЛЯРНЫХ МАРОК ПО ЗАДАННОЙ ПОЛУВАРИОГРАММЕ
!
!############################################################
!     SUBROUTINE VARGEN
!############################################################
!
!
   SUBROUTINE VARGEN(N, NLAG, STP, X, Y, Z, XM, YM, NPR, X1, X10,&
      X11, PAIR, EPS, IS, CPROB, CTEMP, NREFER, NBUC)
!############################################################
      INTEGER N,NLAG,STP,NPR,IS
      INTEGER, DIMENSION(:) :: NBUC
      INTEGER, DIMENSION(:,:) :: NREFER
      REAL*8, DIMENSION(:) :: X(N),Y(N),Z(N), PAIR
      REAL*8, DIMENSION(:,:) :: X1(2,NLAG), X10(2,NLAG), X11(2,NLAG)
      REAL*8 XM,YM,EPS,A,CPROB,CTEMP,RASS,XSTEP,YSTEP,AMEAN,VAR,&
         STDEV,VARCOF,RASS0,RASSK,CPR,TEMP,TEMP0,RATEMP,&
         XN1,YN1,ZN1,XN2,YN2,ZN2,DIS,D1,D2,DRASS,PROB
!============================================================
!----- ВЫЧИСЛЕНИЕ ВАРИОГРАММЫ НАЧАЛЬНОЙ КОНФИГУРАЦИИ --------
      CALL VARIOG(N,NLAG,STP,X,Y,Z,PAIR,X10,AMEAN,VAR,STDEV,VARCOF)
!----- ВЫЧИСЛИМ НАЧАЛЬНУЮ РАЗНОСТЬ --------------------------
      CALL RASVAR(X1,X10,NLAG,STP,RASS0)
      RASSK=RASS0
!----- УПОРЯДОЧЕНИЕ ТОЧЕК -----------------------------------
      RMAX=STP*NLAG
      NX=INT(XM/RMAX)
      NY=INT(YM/RMAX)
      NXY=NX*NY
      CALL BUCK(X,Y,N,XM,YM,NREFER,NBUC,NX,NY,NXY)
      XSTEP=XM/NX
      YSTEP=YM/NY
      CPR=1./CPROB
!----- ПРОЦЕДУРА МЕТРОПОЛИСА --------------------------------
      DO I100=1,21
         IF(I100.EQ.1) THEN
            NPRO=200
         ELSE
            NPRO=NPR
            CPR=CPR*CPROB
            IF(I100.EQ.2) THEN
               TEMP0=SQRT(TEMP/200.)
               TEMP=TEMP0*CTEMP/CPROB
               PRINT 77,TEMP0
77             FORMAT(' НАЧАЛЬНАЯ ТЕМПЕРАТУРА=',F10.5)
            END IF
            TEMP=TEMP*CPROB
            RATEMP=TEMP/TEMP0
            II=I100-1
            PRINT 88,II,RASS0,RATEMP
88          FORMAT(' ШАГ НОМЕР ',I5,' RASS0=',F10.5,&
               ' T/T0=',F10.5)
         END IF
         DO I200=1,NPRO
            DO I70=1,NLAG
               X11(2,I70)=X10(2,I70)
            END DO
!----- ВЫБОР НОМЕРОВ ТОЧЕК ----------------------------------
            CALL GENER(IS, A)
            NOM1=INT(A*N)+1
            CALL GENER(IS, A)
            NOM2=INT(A*N)+1
            XN1=X(NOM1)
            YN1=Y(NOM1)
            ZN1=Z(NOM1)
            XN2=X(NOM2)
            YN2=Y(NOM2)
            ZN2=Z(NOM2)
!----- ОБРАБОТКА ТОЧКИ NOM1 ---------------------------------
!----- ОПРЕДЕЛЕНИЕ НОМЕРА ЦЕНТРАЛЬНОГО BUCKET'А -------------
            IXX=INT(XN1/XSTEP)+1
            IYY=INT(YN1/YSTEP)+1
            DO JY=1,3
               IYC=IYY-2+JY
               IF(IYC.EQ.0.OR.IYC.EQ.NY+1) CYCLE
               DO JX=1,3
                  IXC=IXX-2+JX
                  IF(IXC.EQ.0.OR.IXC.EQ.NY+1) CYCLE
!----- ВЫБОР НОМЕРОВ ТОЧЕК ----------------------------------
                  NOMXY=(IYC-1)*NX+IXC
                  IF(NREFER(2,NOMXY).EQ.0) CYCLE
                  INP=NREFER(1,NOMXY)
                  IENP=INP+NREFER(2,NOMXY)-1
                  DO JT=INP,IENP
                     JJ=NBUC(JT)
                     IF(JJ.EQ.NOM1.OR.JJ.EQ.NOM2) CYCLE
                     DIS=SQRT((XN1-X(JJ))**2+(YN1-Y(JJ))**2)
                     IF(DIS.GE.RMAX) CYCLE
                     NL=INT(DIS/STP)+1
                     D1=(ZN1-Z(JJ))**2*0.5/PAIR(NL)
                     D2=(ZN2-Z(JJ))**2*0.5/PAIR(NL)
                     X11(2,NL)=X11(2,NL)-D1+D2
                  END DO
               END DO
            END DO
!----- ОБРАБОТКА ТОЧКИ NOM2 ---------------------------------
!----- ОПРЕДЕЛЕНИЕ НОМЕРА ЦЕНТРАЛЬНОГО BUCKET'А -------------
            IXX=INT(XN2/XSTEP)+1
            IYY=INT(YN2/YSTEP)+1
            DO JY=1,3
               IYC=IYY-2+JY
               IF(IYC.EQ.0.OR.IYC.EQ.NY+1) CYCLE
               DO JX=1,3
                  IXC=IXX-2+JX
                  IF(IXC.EQ.0.OR.IXC.EQ.NY+1) CYCLE
!----- ВЫБОР НОМЕРОВ ТОЧЕК ----------------------------------
                  NOMXY=(IYC-1)*NX+IXC
                  IF(NREFER(2,NOMXY).EQ.0) CYCLE
                  INP=NREFER(1,NOMXY)
                  IENP=INP+NREFER(2,NOMXY)-1
                  DO JT=INP,IENP
                     JJ=NBUC(JT)
                     IF(JJ.EQ.NOM1.OR.JJ.EQ.NOM2) CYCLE
                     DIS=SQRT((XN2-X(JJ))**2+(YN2-Y(JJ))**2)
                     IF(DIS.GE.RMAX) CYCLE
                     NL=INT(DIS/STP)+1
                     D1=(ZN2-Z(JJ))**2*0.5/PAIR(NL)
                     D2=(ZN1-Z(JJ))**2*0.5/PAIR(NL)
                     X11(2,NL)=X11(2,NL)-D1+D2
                  END DO
               END DO
            END DO
!----- ЗАКОНЧЕНО ПРОБНОЕ ПЕРЕМЕЩЕНИЕ ------------------------
!----- ВЫЧИСЛЕНИЕ НАЧАЛЬНОЙ ТЕМПЕРАТУРЫ ---------------------
            IF(I100.NE.1) THEN
               CALL RASVAR(X1,X11,NLAG,STP,RASS)
               IF(RASS.GT.RASS0) GO TO 400
            ELSE
               CALL RASVAR(X1,X11,NLAG,STP,RASS)
               DRASS=RASS-RASS0
               TEMP=TEMP+DRASS*DRASS
               RASS0=RASS
            END IF
250         CONTINUE
            DO I1=1,NLAG
               X10(2,I)=X11(2,I)
            END DO
            Z(NOM2)=ZN1
            Z(NOM1)=ZN2
            RASS0=RASS
            IF(RASS0.LT.EPS) RETURN
            CYCLE
400         CALL GENER(IS, A)
            PROB=EXP((RASS0-RASS)/TEMP)
            IF (PRO.LT.PROB) GO TO 250
         END DO
      END DO
      RETURN
   END SUBROUTINE VARGEN
!
!
!############################################################
!     SUBROUTINE RASVAR
!############################################################
!
!
   SUBROUTINE RASVAR(X1, X2, NLAG, STP, RASS)
!############################################################
      INTEGER NLAG,STP
      REAL*8, DIMENSION(:,:) :: X1(2,NLAG), X2(2,NLAG)
      REAL*8 RASS
      RASS=0.
      DO I=1,NLAG
         RASS=RASS+ABS(X1(2,I)-X2(2,I))
      END DO
      RASS=RSAS/NLAG
      RETURN
   END SUBROUTINE RASVAR
!
!
!############################################################
!     SUBROUTINE VARIOG
!############################################################
!
!
   SUBROUTINE VARIOG(N, NLAG, STP, X, Y, Z, PAIR, X1,&
      AMEAN, VAR, STDEV, VARCOF)
!############################################################
!     ПРОГРАММА ВЫЧИСЛЕНИЯ ПОЛУВАРИОГРАММЫ
!            (SYDNEY VIERA)
!############################################################
      INTEGER N,NLAG,STP
      REAL*8, DIMENSION(:) :: X(N),Y(N),Z(N),PAIR
      REAL*8, DIMENSION(:,:) :: X1(2,NLAG)
      REAL*8 AMEAN,VAR,STDEV,VARCOF,SUMZ,SUMZ2,XT,YT,ZT,&
         DX,DY,DZ,DHS,DH
      SUMZ=0.0
      SUMZ2=0.0
      DO J=1,NLAG
         PAIR(J)=0.0
         DO I=1,2
            X1(I,J)=0.0
         END DO
      END DO
      DO J=1,N
         SUMZ=SUMZ+Z(J)
         SUMZ2=SUMZ2+Z(J)*Z(J)
      END DO
      AMEAN=SUMZ/FLOAT(N)
      VAR=SUMZ2/FLOAT(N)-AMEAN*AMEAN
      STDEV=SQRT(VAR)
      VARCOF=STDEV/AMEAN
      N1=N-1
      RMAX=FLOAT(NLAG)*STP
      DO K3=1,N1
         XT=X(K3)
         YT=Y(K3)
         ZT=Z(K3)
         L=K3+1
         DO K4=L,N
            DX=XT-X(K4)
            DY=YT-Y(K4)
            DZ=ZT-Z(K4)
            DHS=DX*DX+DY*DY
            IF(DHS.LT.1.0E-6) CYCLE
            DH=SQRT(DHS)
            IF(DH.GT.RMAX) CYCLE
            L4=INT(DH/STP)+1
            PAIR(L4)=PAIR(L4)+1.0
            X1(1,L4)=X1(1,L4)+DH
            X1(2,L4)=X1(2,L4)+0.5*DZ*DZ
         END DO
      END DO
      DO K5=1,NLAG
         IF(PAIR(K5).EQ.0) CYCLE
         X1(1,K5)=X1(1,K5)/PAIR(K5)
         X1(2,K5)=X1(2,L4)/PAIR(K5)
      END DO
      RETURN
   END SUBROUTINE VARIOG
! АЛГОРИТМ ГЕНЕРИРОВАНИЯ РАСШИРЕНИЯ ПРОБНОЙ
! ПЛОЩАДИ С СОХРАНЕНИЕМ ПОЛУВАРИОГРАММЫ
!
!############################################################
!     SUBROUTINE VARGE9
!############################################################
!
!
   SUBROUTINE VARGE9(N, NLAG, STP, X, Y, Z, XM, YM, NPR, X1, X10,&
      X11, PAIR, AMEAN, VAR, STDEV, VARCOF,&
      EPS, IS, CPROB, NREFER, NBUC)
!############################################################
      INTEGER N,NLAG,STP,NPR,IS,NG
      INTEGER, DIMENSION(:) :: NBUC
      INTEGER, DIMENSION(:,:) :: NREFER
      REAL*8, DIMENSION(:) :: X(N),Y(N),Z(N), PAIR
      REAL*8, DIMENSION(:,:) :: X1(2,NLAG), X10(2,NLAG), X11(2,NLAG)
      REAL*8 XM,YM,EPS,A,CPROB,RASS,XSTEP,YSTEP,AMEAN,VAR,&
         STDEV,VARCOF,RASS0,RASSK,CPR,TEMP,TEMP0,RATEMP,&
         XN1,YN1,ZN1,XN2,YN2,ZN2,DIS,D1,D2,DRASS,PROB,CTEMP,&
         XMM,YMM
!============================================================
      DATA CTEMP,CNPR/0.05,3.00/
!----- ВЫЧИСЛЕНИЕ ВАРИОГРАММЫ НАЧАЛЬНОЙ КОНФИГУРАЦИИ --------
      NG=N*9
      XMM=XM*3
      YMM=YM*3
      CALL VARIOG(NG,NLAG,STP,X,Y,Z,PAIR,X10,AMEAN,VAR,STDEV,VARCOF)
!----- ВЫЧИСЛИМ НАЧАЛЬНУЮ РАЗНОСТЬ --------------------------
      CALL RASVAR(X1,X10,NLAG,STP,RASS0)
      RASSK=RASS0
!----- УПОРЯДОЧЕНИЕ ТОЧЕК -----------------------------------
      RMAX=STP*NLAG
      NX=INT(XMM/RMAX)
      NY=INT(YMM/RMAX)
      NXY=NX*NY
      CALL BUCK(X,Y,N,XMM,YMM,NREFER,NBUC,NX,NY,NXY)
      XSTEP=XMM/NX
      YSTEP=YMM/NY
      CPR=1./CPROB
      NPR=INT(NG*CNPR)
!----- ПРОЦЕДУРА МЕТРОПОЛИСА --------------------------------
      DO I100=1,31
         IF(I100.EQ.1) THEN
            NPRO=200
         ELSE
            NPRO=NPR
            CPR=CPR*CPROB
            IF(I100.EQ.2) THEN
               TEMP0=SQRT(TEMP/200.)
               TEMP=TEMP0*CTEMP/CPROB
               PRINT 77,TEMP0
77             FORMAT(' НАЧАЛЬНАЯ ТЕМПЕРАТУРА=',F10.5)
            END IF
            TEMP=TEMP*CPROB
            RATEMP=TEMP/TEMP0
            II=I100-1
            PRINT 88,II,RASS0,RATEMP
88          FORMAT(' ШАГ НОМЕР ',I5,' RASS0=',F10.5,&
               ' T/T0=',F10.5)
         END IF
         DO I200=1,NPRO
            DO I70=1,NLAG
               X11(2,I70)=X10(2,I70)
            END DO
!----- ВЫБОР НОМЕРОВ ТОЧЕК ----------------------------------
            CALL GENER(IS, A)
            NOM1=N+INT(A*(NG-N))+1
            CALL GENER(IS, A)
            NOM2=N+INT(A*(NG-N))+1
            XN1=X(NOM1)
            YN1=Y(NOM1)
            ZN1=Z(NOM1)
            XN2=X(NOM2)
            YN2=Y(NOM2)
            ZN2=Z(NOM2)
!----- ОБРАБОТКА ТОЧКИ NOM1 ---------------------------------
!----- ОПРЕДЕЛЕНИЕ НОМЕРА ЦЕНТРАЛЬНОГО BUCKET'А -------------
            IXX=INT(XN1/XSTEP)+1
            IYY=INT(YN1/YSTEP)+1
            DO JY=1,3
               IYC=IYY-2+JY
               IF(IYC.EQ.0.OR.IYC.EQ.NY+1) CYCLE
               DO JX=1,3
                  IXC=IXX-2+JX
                  IF(IXC.EQ.0.OR.IXC.EQ.NY+1) CYCLE
!----- ВЫБОР НОМЕРОВ ТОЧЕК ----------------------------------
                  NOMXY=(IYC-1)*NX+IXC
                  IF(NREFER(2,NOMXY).EQ.0) CYCLE
                  INP=NREFER(1,NOMXY)
                  IENP=INP+NREFER(2,NOMXY)-1
                  DO JT=INP,IENP
                     JJ=NBUC(JT)
                     IF(JJ.EQ.NOM1.OR.JJ.EQ.NOM2) CYCLE
                     DIS=SQRT((XN1-X(JJ))**2+(YN1-Y(JJ))**2)
                     IF(DIS.GE.RMAX) CYCLE
                     NL=INT(DIS/STP)+1
                     D1=(ZN1-Z(JJ))**2*0.5/PAIR(NL)
                     D2=(ZN2-Z(JJ))**2*0.5/PAIR(NL)
                     X11(2,NL)=X11(2,NL)-D1+D2
                  END DO
               END DO
            END DO
!----- ОБРАБОТКА ТОЧКИ NOM2 ---------------------------------
!----- ОПРЕДЕЛЕНИЕ НОМЕРА ЦЕНТРАЛЬНОГО BUCKET'А -------------
            IXX=INT(XN2/XSTEP)+1
            IYY=INT(YN2/YSTEP)+1
            DO JY=1,3
               IYC=IYY-2+JY
               IF(IYC.EQ.0.OR.IYC.EQ.NY+1) CYCLE
               DO JX=1,3
                  IXC=IXX-2+JX
                  IF(IXC.EQ.0.OR.IXC.EQ.NY+1) CYCLE
!----- ВЫБОР НОМЕРОВ ТОЧЕК ----------------------------------
                  NOMXY=(IYC-1)*NX+IXC
                  IF(NREFER(2,NOMXY).EQ.0) CYCLE
                  INP=NREFER(1,NOMXY)
                  IENP=INP+NREFER(2,NOMXY)-1
                  DO JT=INP,IENP
                     JJ=NBUC(JT)
                     IF(JJ.EQ.NOM1.OR.JJ.EQ.NOM2) CYCLE
                     DIS=SQRT((XN2-X(JJ))**2+(YN2-Y(JJ))**2)
                     IF(DIS.GE.RMAX) CYCLE
                     NL=INT(DIS/STP)+1
                     D1=(ZN2-Z(JJ))**2*0.5/PAIR(NL)
                     D2=(ZN1-Z(JJ))**2*0.5/PAIR(NL)
                     X11(2,NL)=X11(2,NL)-D1+D2
                  END DO
               END DO
            END DO
!----- ЗАКОНЧЕНО ПРОБНОЕ ПЕРЕМЕЩЕНИЕ ------------------------
!----- ВЫЧИСЛЕНИЕ НАЧАЛЬНОЙ ТЕМПЕРАТУРЫ ---------------------
            IF(I100.NE.1) THEN
               CALL RASVAR(X1,X11,NLAG,STP,RASS)
               IF(RASS.GT.RASS0) GO TO 400
            ELSE
               CALL RASVAR(X1,X11,NLAG,STP,RASS)
               DRASS=RASS-RASS0
               TEMP=TEMP+DRASS*DRASS
               RASS0=RASS
            END IF
250         CONTINUE
            DO I1=1,NLAG
               X10(2,I)=X11(2,I)
            END DO
            Z(NOM2)=ZN1
            Z(NOM1)=ZN2
            RASS0=RASS
            IF(RASS0.LT.EPS) GO TO 2000
            CYCLE
400         CALL GENER(IS, A)
            PROB=EXP((RASS0-RASS)/TEMP)
            IF (PRO.LT.PROB) GO TO 250
         END DO
      END DO
2000  N=NG
      XM=XMM
      YM=YMM
      PRINT 741,RASS0
741   FORMAT(' RASS0=',F12.6)
      RETURN
   END SUBROUTINE VARGE9


end module GEN

