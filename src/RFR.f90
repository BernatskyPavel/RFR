module RFR
   private

   PUBLIC :: GENRFR
contains
   
! Алгоритм генерирования размещений по радиальной функции распределения
!
!############################################################
!     SUBROUTINE
!############################################################
!
!
   SUBROUTINE GENRFR(N, X, Y, XM, YM, RMAX, K1,&
      RF, RK, S, RF0, S0, RF1, S1, EPSR,&
      NBUC, NREFER, IS, NPR, CPROB)
!############################################################
      REAL*8, DIMENSION (:) :: X(N),Y(N),RF(K1),RK(K1),S(K1),&
         RF1(K1),S1(K1),RF0(K1),S0(K1),&
         RF0N(200),RF1N(200)
      INTEGER, DIMENSION(:) :: NBUC
      INTEGER, DIMENSION(:,:) :: NREFER
      REAL*8 A,PRO,PROB,RASS0,RASSK,RH2,RR2,ORH2,&
         XM,YM,RMAX,DENSE,XSTEP,YSTEP,DRASS,DIS2,&
         XNOM,YNOM,TEMP,RAD,TEMP0,RATEMP,XJ,YJ,RASS,&
         EPSR,CPROB,CPR
!============================================================
      DATA CRMA,CTEMP/4.,0.5/
!------- ГЕНЕРИРОВАНИЕ НАЧАЛЬНОГО РАЗМЕЩЕНИЯ ----------------
      CALL SRAND1()
      DO I=1,N
         CALL GENER(IS,A)
         X(I)=A*XM
         CALL GENER(IS,A)
         Y(I)=A*YM
      END DO

      OPEN(UNIT=48,FILE='data_x.txt')
      DO I=1,N
         WRITE(48,*) X(I),Y(I)
      END DO
      CLOSE(48)

      CALL SYSTEM('gnuplot -p data_plot_x.plt')

      RH2=RK(1)*RK(1)
      K2=INT(RMAX*RMAX/RH2)
      IF(K2.GT.K1) GO TO 15
      K1=K2
15    RMAX=SQRT(RH2*K1)
      RR2=RH2*K1
      ORH2=1./RH2
!------ ВЫЧИСЛЕНИЕ RFR НАЧАЛЬНОГО РАЗМЕЩЕНИЯ ----------------
      ALF=0.1
      EPS=1.96
      CALL RFR1(N,XM,YM,X,Y,ALF,EPS,RMAX,RF0,RK,S0,K1,1)
      DENSE=FLOAT(N)/XM/YM
      DO I=1,K1
         RF0N(I) = RF0(I)*S0(I)*DENSE
         RK(I)=RH2*I
      END DO
!------ ВЫЧИСЛИМ НАЧАЛЬНУЮ РАЗНОСТЬ -------------------------
      CALL RASST(RF,RF0,RK,K1,RASS0)
      RASSK=RASS0
!------------- УПОРЯДОЧЕНИЕ ТОЧЕК ---------------------------
      NX=INT(XM/RMAX)
      NY=INT(YM/RMAX)
      NXY=NX*NY
      CALL BUCK(X,Y,N,XM,YM,NREFER,NBUC,NX,NY,NXY)
      XSTEP=XM/FLOAT(NX)
      YSTEP=YM/FLOAT(NY)
      TEMP=0.
      CPR=1./CPROB
!------ ПРОЦЕДУРА МЕТРОПОЛИСА -------------------------------
      DO ITP=1,61
         RAD=CRMA*RMAX
         IF (ITP.NE.1) GO TO 17
         NPRO=200
         GO TO 18
17       NPRO=NPR
         CPR=CPR*CPROB
         IF(ITP.EQ.2) GO TO 171
         TEMP0=SQRT(TEMP/200.)
         TEMP=TEMP0*CTEMP/CPROB
         PRINT 77,TEMP0
77       FORMAT(' НАЧАЛЬНАЯ ТЕМПЕРАТУРА =',F10.5)
171      TEMP=TEMP*CPROB
         RATEMP=TEMP/TEMP0
         III=ITP-1
         PRINT 88, III, RASS0, RATEMP
88       FORMAT(' ШАГ НОМЕР ',I5,'  RASS0=', F10.5, 'T/T0=', F10.6)
18       CONTINUE
         DO ITP1=1,NPRO
            DO I=1,K1
               RF1(I)=RF0(I)
               RF1N(I)=RF0N(I)
               S1(I)=S0(I)
            END DO
            !---------- ВЫБОР НОМЕРА ТОЧКИ ------------------------------
            CALL GENER(IS,A)
            NOM=INT(A*N)+1
            !--------- ОПРЕДЕЛЕНИЕ СОСЕДЕЙ И РАССТОЯНИЯ ДО НИХ ----------
            XJ=X(NOM)
            YJ=Y(NOM)
            !--------- ОПРЕДЕЛЕНИЕ НОМЕРА ЦЕНТРАЛЬНОГО BUCKET'А ---------
            IXX=INT(XJ/XSTEP)+1
            IYY=INT(YJ/YSTEP)+1
            DO JY=1,3
               IYC=IYY-2+JY
               IF(IYC.EQ.0.OR.IYC.EQ.NY+1) CYCLE
               DO JX=1,3
                  IXC=IXX-2+JX
                  IF(IXC.EQ.0.OR.IXC.EQ.NX+1) CYCLE
!---------- ВЫБОР НОМЕРОВ ТОЧЕК -----------------------------
                  NOMXY=(IYC-1)*NX+IXC
                  IF(NREFER(2,NOMXY).EQ.0) CYCLE
                  INP=NREFER(1,NOMXY)
                  IENP=INP+NREFER(2,NOMXY)-1
                  DO JT=INP,IENP
                     JJ=NBUC(JT)!BUG =(
                     IF(JJ.EQ.NOM) CYCLE
                     !IF(JJ.GT.0)
                     DIS2=(XJ-X(JJ))**2+(YJ-Y(JJ))**2
                     IF(DIS2.GE.RR2) CYCLE
                     NL=INT(DIS2*ORH2)+1
                     RF1N(NL)=RF1N(NL)-2
                  END DO
               END DO
            END DO

            CALL AREA1(XM,YM,XJ,YJ,K1,RK,S1,RH2,1)
!---------- ТОЧКА ИЗ`ЯТА, RF1 и S1 СКОРРЕКТИРОВАНЫ ----------
!---------- ПОПРОБУЕМ ПЕРЕМЕСТИТЬ ТОЧКУ NOM -----------------
            CALL GENER(IS,A)
            XNOM=X(NOM)+(A-0.5)*RAD
            XNOM=XNOM-FLOAT(IENT(XNOM/XM))*XM
            CALL GENER(IS,A)
            YNOM=Y(NOM)+(A-0.5)*RAD
            YNOM=YNOM-FLOAT(IENT(YNOM/YM))*YM
!--------- ОПРЕДЕЛЕНИЕ НОМЕРА ЦЕНТРАЛЬНОГО BUCKET'А ---------
            IXX1=INT(XNOM/XSTEP)+1
            IYY1=INT(YNOM/YSTEP)+1
            DO JY=1,3
               IYC=IYY1-2+JY
               IF(IYC.EQ.0.OR.IYC.EQ.NY+1) CYCLE
               DO JX=1,3
                  IXC=IXX1-2+JX
                  IF(IXC.EQ.0.OR.IXC.EQ.NX+1) CYCLE
!---------- ВЫБОР НОМЕРОВ ТОЧЕК -----------------------------
                  NOMXY1=(IYC-1)*NX+IXC
                  IF(NREFER(2,NOMXY1).EQ.0) CYCLE
                  INP=NREFER(1,NOMXY1)
                  IENP=INP+NREFER(2,NOMXY1)-1
                  !IF(INP.LE.0) INP=1
                  !IF(IENP.GT.N) IENP=N
                  DO JT=INP,IENP
                     JJ=NBUC(JT)
                     IF(JJ.EQ.NOM) CYCLE
                     DIS2=(XNOM-X(JJ))**2+(YNOM-Y(JJ))**2
                     IF(DIS2.GE.RR2) CYCLE
                     NL=INT(DIS2*ORH2)+1
                     RF1N(NL)=RF1N(NL)+2.
                  END DO
               END DO
            END DO
            CALL AREA1(XM,YM,XNOM,YNOM,K1,RK,S1,RH2,0)
!--------- ЗАКОНЧЕНО ПРОБНОЕ ПЕРЕМЕЩЕНИЕ --------------------
            DO L=1,K1
               RF0(L)=RF0N(L)/(S0(L)*DENSE)
               RF1(L)=RF1N(L)/(S1(L)*DENSE)
            END DO
!--------- ВЫЧИСЛЕНИЕ НАЧАЛЬНОЙ ТЕМПЕРАТУРЫ -----------------
            IF(ITP.NE.1) GO TO 220
            CALL RASST(RF,RF1,RK,K1,RASS)
            DRASS=RASS-RASS0
            TEMP=TEMP+DRASS*DRASS
            RASS0=RASS
            GO TO 250
220         CALL RASST(RF,RF1,RK,K1,RASS)
            IF(RASS.GT.RASS0) GO TO 400
!--------- РАЗНИЦА УМЕНЬШИЛАСЬ ! ----------------------------
250         CONTINUE
            DO I=1,K1
               RF0N(I)=RF1N(I)
               S0(I)=S1(I)
            END DO
            X(NOM)=XNOM
            Y(NOM)=YNOM
            RASS0=RASS
!--------- ПЕРЕУПОРЯДОЧИВАНИЕ ПОСЛЕ УДАЧНОГО СДВИГА ---------
            NOLD=(IYY-1)*NX+IXX
            NNEW=(IYY1-1)*NX+IXX1
            I1=NREFER(1,NOLD)
            !IF(I1.LE.0) I1=1
            I2=I1+NREFER(2,NOLD)-1
            !IF(I2.GT.N) I2=N
            DO L=I1,I2
               IF(NBUC(L).NE.NOM) CYCLE
               NOMB1=L
               EXIT
            END DO
            NREFER(2,NOLD)=NREFER(2,NOLD)-1
            NREFER(2,NNEW)=NREFER(2,NNEW)+1
            IF(NNEW.LE.NOLD) GO TO 600
            NOMB2=NREFER(1,NNEW)-1
            !IF(NOMB.LT.1) NOMB=1
            DO L=NOMB1,NOMB2
               IF(L.LT.N) NBUC(L)=NBUC(L+1)
            END DO
            !IF(NOMB2.GT.N) NOMB2=N
            !IF(NOMB2.GT.0)
            NBUC(NOMB2)=NOM
            NU=NOLD+1
            !IF(NU.LE.0) NU = 1
            !IF(NNEW.GT.NXY) NNEW=NXY
            DO L=NU,NNEW
               NREFER(1,L)=NREFER(1,L)-1
            END DO
            CYCLE
600         IF(NNEW.EQ.NOLD) GO TO 700
            NOMB2=NREFER(1,NNEW)
            LZ=NOMB1-NOMB2
            DO L=1,LZ
               IL=NOMB1-L+1
               !IF(IL.GT.1.AND.IL.LE.N)
               NBUC(IL)=NBUC(IL-1)
            END DO
            !IF(NOMB2.GT.N) NOMB2=N
            !IF(NOMB2.GT.0)
            NBUC(NOMB2)=NOM
            NU=NNEW+1
            DO L=NU,NOLD
               NREFER(1,L)=NREFER(1,L)+1
            END DO
700         IF(RASS0.LT.EPSR) GO TO 2000
            CYCLE
400         CALL GENER(IS,PRO)
            PROB=EXP((RASS0-RASS)/TEMP)
            IF(PRO.LT.PROB) GO TO 250
         END DO
      END DO
      !------- КОНЕЦ ЦИКЛА ПРОЦЕДУРЫ МЕТРОПОЛИСА ------------------
2000  CONTINUE
      DO I=1,K1
         RK(I) = SQRT(RK(I))
      END DO
      RETURN
   END SUBROUTINE
   !
   !
   !############################################################
   !     SUBROUTINE BUCK
   !############################################################
   !
   !
   SUBROUTINE BUCK(X,Y,N,XM,YM,NREFER,NBUC,NX,NY,NXY)
      REAL*8, DIMENSION(:) :: X(N),Y(N)
      INTEGER NX,NY,NXY,NBUC(N),NREFER(2,NXY)
      REAL*8 XM,YM,DX,DY
      DX=XM/FLOAT(NX)
      DY=YM/FLOAT(NY)
      !------- ПЕРВОЕ ЗАПОЛНЕНИЕ МАССИВА NREFFER(2, ) -------------
      DO I=1,N
         IX=INT(X(I)/DX)+1
         IY=INT(Y(I)/DY)+1
         J=(IY-1)*NX+IX
         NREFER(2,J)=NREFER(2,J)+1
      END DO
      !------- ЗАПОЛНЕНИЕ ПЕРВОЙ СТРОКИ NREFFER - НОМЕР В NBUC ----
      ITIC=1
      DO I=1,NXY
         NREFER(1,I)=ITIC
         ITIC=NREFER(2,I)+ITIC
      END DO
      DO I=1,NXY
         NREFER(2,I)=0
      END DO
      !------- ЗАПОЛНЕНИЕ NBUC И 2-ОЙ СТРОКИ NREFFER --------------
      DO I=1,N
         IX=INT(X(I)/DX)+1
         IY=INT(Y(I)/DY)+1
         NOM=(IY-1)*NX+IX
         NOMXY=NREFER(1,NOM)+NREFER(2,NOM)
         !IF(NOMXY.LE.N)
         NBUC(NOMXY)=I
         NREFER(2,NOM)=NREFER(2,NOM)+1
      END DO
      RETURN
   END SUBROUTINE
   !
   !
   !############################################################
   !     SUBROUTINE RFR
   !############################################################
   !
   !
   SUBROUTINE RFR1(N1,XX,YY,X,Y,ALF,EPS,RMAX,RF,RK,&
      S,K1,KEY)
      !############################################################
      REAL*8, DIMENSION(:) :: X(N1),Y(N1),RF(K1),RK(K1),S(K1)
      REAL*8 XX,YY,PLS,RMAX,RH2,RS,XX2,YY2,RR2,RR,ORH2,R2,XI,YI
      PLS=N1/(XX*YY)
      RS=SQRT(1./PLS)
      RH2=RK(1)*RK(1)
      IF(KEY.EQ.1) GO TO 5
      H1=EPS/(ALF*SQRT(2.*ASIN(1.)*N1))
      RH2=H1*H1/PLS
5     XX2=XX/2
      YY2=YY/2
      IF(KEY.EQ.1) GO TO 15
      K1=INT(((AMIN1(XX2,YY2,RMAX))**2)/RH2)
      PRINT 201,N1,XX,YY,EPS,ALF
201   FORMAT (//1X,'ТОЧЕК',I5,'УЧАСТОК',2F5.0,'EPS=',F7.2,&
         'ALF=',F7.4)
      PRINT 200,RS,PLS
200   FORMAT (1X,'СРЕДНЕЕ РАССТОЯНИЕ',F7.3,&
         'ПЛОТНОСТЬ',F7.3)
      PRINT 202,H1,K1
202   FORMAT (1X,'ШАГ',F7.5,'ЧИСЛО ШАГОВ',I5)
15    ORH2=1./RH2
      RR2=RH2*K1
      RR=SQRT(RR2)
      DO I=1,K1
         RK(I)=RH2*I
         RF(I)=0.
      END DO
      !     CALL SHELL(N1,X,Y)
      I2=N1-1
      DO I=1,I2
         XI=X(I)
         YI=Y(I)
         J1=I+1
         DO J=J1,N1
            IF((X(J)-X(I)).GE.RR) EXIT
            R2=(X(J)-XI)**2+(Y(J)-YI)**2
            IF(R2.GE.RR2) CYCLE
            N=INT(R2*ORH2)+1
            RF(N)=RF(N)+2.
         END DO
      END DO
      CALL AREA(N1,XX,YY,X,Y,K1,RK,S,RH2)
      DO I=1,K1
         RF(I)=RF(I)/(S(I)*PLS)
         RK(I)=SQRT(RK(I))
      END DO
      PRINT 205
205   FORMAT (1X,'РАДИАЛЬНАЯ ФУНКЦИЯ'&
         'РАСПРЕДЕЛЕНИЯ')
      !     CALL DRAW(K1,RK,RF)
      J1=1
      J2=10
      DO I=1,K1,10
         IF(J2.GT.K1) J2=K1
         PRINT 224,(RF(J),J=J1,J2)
         J1=J1+10
         J2=J2+10
      END DO
224   FORMAT (1X,'RFR',10F11.4)
      RETURN
   END SUBROUTINE
   !
   !
   !############################################################
   !     SUBROUTINE AREA1
   !############################################################
   !
   !
   SUBROUTINE AREA1(XX,YY,X,Y,K1,RK,S,RH2,KEY)
      REAL*8, DIMENSION(:) :: RK(K1),S(K1)
      REAL*8 S01,XX,YY,H,G,X,Y,RH2,ORH2,RR2,RR
      ORH2=1./RH2
      RR2=RH2*K1
      RR=SQRT(RR2)
      S01=RK(1)*2.*ASIN(1.)
      IF(KEY.EQ.0) GO TO 2
      DO I=1,K1
         S(I)=S(I)-S01
      END DO
      GO TO 6
2     CONTINUE
      DO I=1,K1
         S(I)=S(I)+S01
      END DO
6     CONTINUE
      IF(X.GE.RR) GO TO 20
      CALL SEGM1(K1,RK,S,ORH2,X,KEY)
      IF(Y.GE.RR) GO TO 11
      CALL SEGM1(K1,RK,S,ORH2,Y,KEY)
      IF((X*X+Y*Y).GE.R2) RETURN
      CALL CROSS1(K1,RK,S,ORH2,X,Y,KEY)
      RETURN
11    G=YY-Y
      IF(G.GE.RR) RETURN
      CALL SEGM1(K1,RK,S,ORH2,G,KEY)
      IF((X*X+G*G).GE.RR2) RETURN
      CALL CROSS1(K1,RK,S,ORH2,X,G,KEY)
      RETURN
20    H=XX-X
      IF(H.LT.RR) GO TO 30
      IF(Y.GE.RR) GO TO 21
      CALL SEGM1(K1,RK,S,ORH2,Y,KEY)
      RETURN
21    G=YY-Y
      IF(G.GE.RR) RETURN
      CALL SEGM1(K1,RK,S,ORH2,G,KEY)
      RETURN
30    CALL SEGM1(K1,RK,S,ORH2,H,KEY)
      IF(Y.GE.RR) GO TO 31
      CALL SEGM1(K1,RK,S,ORH2,Y,KEY)
      IF((H*H+Y*Y).GE.RR2) RETURN
      CALL CROSS1(K1,RK,S,ORH2,H,Y,KEY)
31    G=YY-Y
      IF(G.GE.RR) RETURN
      CALL SEGM1(K1,RK,S,ORH2,G,KEY)
      IF((H*H+G*G).GE.RR2) RETURN
      CALL CROSS1(K1,RK,S,ORH2,H,G,KEY)
      RETURN
   END SUBROUTINE
   !
   !
   !############################################################
   !     SUBROUTINE SEGM1
   !############################################################
   !
   !
   SUBROUTINE SEGM1(K1,RK,S,ORH2,X,KEY)
      REAL*8, DIMENSION(:) :: RK(K1),S(K1)
      REAL*8 C,P,ORH2,X,X2
      X2=X*X
      KX=INT(X2*ORH2)+1
      DO I=KX,K1-1
         C=SQRT(RK(I)-X2)
         P=RK(I)*ATAN(C/X)-C*X
         IF(KEY.EQ.1) GO TO 2
         S(I)=S(I)-P
         S(I+1)=S(I+1)+P
         CYCLE
2        S(I)=S(I)+P
         S(I+1)=S(I+1)-P
      END DO
      RETURN
   END SUBROUTINE
   !
   !
   !############################################################
   !     SUBROUTINE CROSS1
   !############################################################
   !
   !
   SUBROUTINE CROSS1(K1,RK,S,ORH2,X,Y,KEY)
      REAL*8, DIMENSION(:) :: RK(K1),S(K1)
      REAL*8 A,B,P,C,X,Y,X2,Y2,ORH2
      X2=X*X
      Y2=Y*Y
      KXY=INT((X2+Y2)*ORH2)+1
      DO I=KXY,K1-1
         A=SQRT(RK(I)-X2)-Y
         B=SQRT(RK(I)-Y2)-X
         C=0.25*(A*A+B*B)
         P=0.5*A*B-RK(I)*ASIN(SQRT(C/RK(I)))-SQRT(C*(RK(I)-C))
         IF(KEY.EQ.1) GO TO 2
         S(I)=S(I)+P
         S(I+1)=s(I+1)-P
         CYCLE
2        S(I)=S(I)-P
         S(I+1)=S(I+1)+P
      END DO
      RETURN
   END SUBROUTINE
   !
   !
   !############################################################
   !     SUBROUTINE RASST
   !############################################################
   !
   !
   SUBROUTINE RASST(RF,RF1,RK,K1,RASS)
      REAL*8, DIMENSION (:) :: RF(K1),RF1(K1),RK(K1)
      REAL*8 W, RASS
      RASS=0.
      DO I=1,K1
         IF(I.GT.1) GO TO 2
         W=SQRT(RK(1))
         GO TO 3
2        W=(SQRT(RK(I))-SQRT(RK(I-1)))
3        RASS=RASS+ABS(RF(I)-RF1(I))*W
      END DO
      RASS=RASS/SQRT(RK(K1))
      RETURN
   END SUBROUTINE
   !
   !
   !############################################################
   !     SUBROUTINE AREA
   !############################################################
   !
   !
   SUBROUTINE AREA(N1,XX,YY,X,Y,K1,RK,S,RH2)
      REAL*8, DIMENSION(:) :: X(N1),Y(N1),RK(K1),S(K1)
      REAL*8 G,H,XX,YY,RH2,ORH2,RR2,RR,S0
      ORH2=1./RH2
      S0=2.*ASIN(1.)*RH2*N1
      DO J=1,K1
         S(J)=S0
      END DO
      RR2=RH2*K1
      RR=SQRT(RR2)
      I=1
10    IF(X(I).GE.RR) GO TO 20
      CALL SEGM(K1,RK,S,ORH2,X(I))
      IF(Y(I).GE.RR) GO TO 11
      CALL SEGM(K1,RK,S,ORH2,Y(I))
      IF((X(I)*X(I)+Y(I)*Y(I)).GE.RR2) GO TO 11
      CALL CROSS(K1,RK,S,ORH2,X(I),Y(I))
      GO TO 12
11    G=YY-Y(I)
      IF(G.GE.RR) GO TO 12
      CALL SEGM(K1,RK,S,ORH2,G)
      IF((X(I)*X(I)+G*G).GE.RR2) GO TO 12
      CALL CROSS(K1,RK,S,ORH2,X(I),G)
12    I=I+1
      IF(I.GT.N1) RETURN
      GO TO 10
20    H=XX-X(I)
      IF(H.LT.RR) GO TO 30
      IF(Y(I).GE.RR) GO TO 21
      CALL SEGM(K1,RK,S,ORH2,Y(I))
      GO TO 22
21    G=YY-Y(I)
      IF(G.GE.RR) GO TO 22
      CALL SEGM(K1,RK,S,ORH2,G)
22    I=I+1
      IF(I.GT.N1) RETURN
      GO TO 20
30    CALL SEGM(K1,RK,S,ORH2,H)
      IF(Y(I).GE.RR) GO TO 31
      CALL SEGM(K1,RK,S,ORH2,Y(I))
      IF((H*H+Y(I)*Y(I)).GE.RR2) GO TO 31
      CALL CROSS(K1,RK,S,ORH2,H,Y(I))
31    G=YY-Y(I)
      IF(G.GE.RR) GO TO 32
      CALL SEGM(K1,RK,S,ORH2,G)
      IF((H*H+G*G).GE.RR2) GO TO 32
      CALL CROSS(K1,RK,S,ORH2,H,G)
32    I=I+1
      IF(I.GT.N1) RETURN
      H=XX-X(I)
      GO TO 30
   END SUBROUTINE
   !
   !
   !############################################################
   !     SUBROUTINE SEGM
   !############################################################
   !
   !
   SUBROUTINE SEGM(K1,RK,S,ORH2,X)
      REAL*8, DIMENSION(:) :: RK(K1),S(K1)
      REAL*8 C,P,X,X2,ORH2
      X2=X*X
      KX=INT(X2*ORH2)+1
      DO I=KX,K1-1
         C=SQRT(RK(I)-X2)
         P=RK(I)*ATAN(C/X)-C*X
         S(I)=S(I)-P
         S(I+1)=S(I+1)+P
      END DO
      RETURN
   END SUBROUTINE
   !
   !
   !############################################################
   !     SUBROUTINE CROSS
   !############################################################
   !
   !
   SUBROUTINE CROSS(K1,RK,S,ORH2,X,Y)
      REAL*8, DIMENSION(:) :: RK(K1),S(K1)
      REAL*8 X,Y,A,B,C,P,X2,Y2,ORH2
      X2=X*X
      Y2=Y*Y
      KXY=INT((X2+Y2)*ORH2)+1
      DO I=KXY,K1-1
         A=SQRT(RK(I)-X2)-Y
         B=SQRT(RK(I)-Y2)-X
         C=0.25*(A*A+B*B)
         P=0.5*A*B-RK(I)*ASIN(SQRT(C/RK(I)))&
            -SQRT(C*(RK(I)-C))
         S(I)=S(I)-P
         S(I+1)=S(I+1)+P
      END DO
      RETURN
   END SUBROUTINE
   !
   !
   !############################################################
   !     SUBROUTINE IENT
   !############################################################
   !
   !
   INTEGER FUNCTION IENT(X)
      REAL*8 X
      IF(X.LT.0) THEN
         IENT=INT(X-1.)
      ELSE
         IENT=INT(X)
      END IF
   END FUNCTION
   !
   !
   !############################################################
   !     SUBROUTINE IENT
   !############################################################
   !
   !
   REAL*8 FUNCTION SIGN1(X)
      IF(X.GT.0) THEN
         SIGN1=1.
      ELSEIF(X.LT.0) THEN
         SIGN1=-1.
      ELSE
         SIGN1=0.
      END IF
      RETURN
   END FUNCTION
   !
   !
   !############################################################
   !     SUBROUTINE GENER
   !############################################################
   !
   !
   SUBROUTINE GENER(IS, A)
      INTEGER IS
      REAL*8 A
      A = RAND()
      !   IS=3125*IS
      !   IS=IS-INT(SIGN1(IS/67108864.)*67108864.)
      !   A=IS/67108864
      RETURN
   END SUBROUTINE

   SUBROUTINE SRAND1()
      INTEGER, ALLOCATABLE :: SEED(:)
      INTEGER :: N
      CALL RANDOM_SEED(SIZE=N)
      ALLOCATE(SEED(N))
      CALL RANDOM_SEED(GET=SEED)
      CALL SRAND(SEED(N))
   END SUBROUTINE


end module RFR
