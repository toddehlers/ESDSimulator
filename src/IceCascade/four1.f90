! FOUR1

! this routine is taken directly out of numerical recipes
! see the book for further information

! subroutines called:
! NONE

        module m_four1
          contains
          SUBROUTINE FOUR1(DATA,NN,ISIGN)
          implicit none
          REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
          real(8) :: DATA, TEMPI, TEMPR

          integer(4) :: NN, ISIGN, I, J, MMAX, N, M, ISTEP

          DIMENSION DATA(*)

          N=2*NN
          J=1
          DO 11 I=1,N,2
            IF(J.GT.I)THEN
              TEMPR=DATA(J)
              TEMPI=DATA(J+1)
              DATA(J)=DATA(I)
              DATA(J+1)=DATA(I+1)
              DATA(I)=TEMPR
              DATA(I+1)=TEMPI
            ENDIF
            M=N/2
    1       IF ((M.GE.2).AND.(J.GT.M)) THEN
              J=J-M
              M=M/2
            GO TO 1
            ENDIF
            J=J+M
11        CONTINUE
          MMAX=2
    2     IF (N.GT.MMAX) THEN
            ISTEP=2*MMAX
            THETA=6.28318530717959D0/(dble(ISIGN*MMAX))
            WPR=-2.D0*DSIN(0.5D0*THETA)**2
            WPI=DSIN(THETA)
            WR=1.D0
            WI=0.D0
            DO 13 M=1,MMAX,2
              DO 12 I=M,N,ISTEP
                J=I+MMAX
                TEMPR=dble(WR)*DATA(J)-dble(WI)*DATA(J+1)
                TEMPI=dble(WR)*DATA(J+1)+dble(WI)*DATA(J)
                DATA(J)=DATA(I)-TEMPR
                DATA(J+1)=DATA(I+1)-TEMPI
                DATA(I)=DATA(I)+TEMPR
                DATA(I+1)=DATA(I+1)+TEMPI
12              CONTINUE
              WTEMP=WR
              WR=WR*WPR-WI*WPI+WR
              WI=WI*WPR+WTEMP*WPI+WI
13            CONTINUE
            MMAX=ISTEP
          GO TO 2
          ENDIF
          RETURN
          END subroutine four1
        end module m_four1


