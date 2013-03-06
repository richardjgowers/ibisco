

! SUBROUTINE PER LA PARSERIZZAZIONE DELLE STRINGHE
!=================================================

        SUBROUTINE PARSE()
	USE MODULEPARSING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!     Parse input string.                                              !
!                                                                      !
!     G. Kneller                                                       !
!     IBM Corp., Data Systems Division, Dept. 48B                      !
!     Kingston, NY, USA                                                !
!                                                                      !
!     ARGUMENTS:                                                       !
!                                                                      !
!     LINE    : Input string to be parsed.                    (INPUT)  !
!               >> character!80 LINE <<                                !
!     SEP     : List with separation characters.              (INPUT)  !
!               >> character!1 SEP(NSEP)                               !
!     NSEP    : Number of separation characters.              (INPUT)  !
!               >> integer NSEP <<                                     !
!     COMM    : Labels for begin and end of comments.         (INPUT)  !
!               >> character!1 COMM(2) <<                              !
!               COMM(1) Label for begin of comment.                    !
!               COMM(2) Label for end of comment.                      !
!     STRNGS  : Extracted strings.                            (OUTPUT) !
!               >> character!80 STRNGS(NSTRNG) <<                      !
!     NSTRNG  : Dimension of STRNGS.                          (INPUT)  !
!               >> integer NSTRNG <<                                   !
!     NWORD   : Number of extracted strings                   (OUTPUT) !
!               >> integer NWORD <<                                    !
!     IRET    : Return code.                                  (OUTPUT) !
!               >> integer IRET <<                                     !
!               IRET = 0,1. For each value of IRET a                   !
!               corresponding error message is assigned to             !
!               ERRMSG.                                                !
!     ERRMSG  : Error message.                                (OUTPUT) !
!               >> character!80 ERRMSG <<                              !
!                                                                      !
!---- LAST UPDATE: 03/07/89 -------------------------------------------!
!                                                                      !
!     Written by Gerald Kneller Dept 48B, IBM Kingston 1989            !
!                                                                      !
!     EXTERNALS: none                                                  !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!==== DECLARATIONS: ===================================================*

      implicit CHARACTER*80(a-z)

!---- PARAMETERS ------------------------------------------------------*

      integer linlen,loclen
      parameter(linlen=80,loclen=80)

!---- ARGUMENTS: ------------------------------------------------------*


!*---- LOCAL VARIABLES: ------------------------------------------------*

      integer      i,j,isep,offset,iword,ifind
      integer      locb(loclen),loce(loclen)
      character*1  linarr(linlen)
      character*80 blank
      logical      sepchr

!---- DATA STATEMENTS: ------------------------------------------------*

      data blank /'                              '/

!==== EXECUTABLE STATEMENTS: ==========================================*

!---- INITIALIZE: -----------------------------------------------------*

      offset = 0
      nword  = 0
      do 100 i = 1,nstrng
         strngs(i) = blank
  100 continue
      do 110 i = 1,linlen
         linarr(i) = ' '
  110 continue

      read(line,'(80a1)',END=1000)(linarr(i),i=1,linlen)
      linarr(80)=' '

!---- PARSE: ----------------------------------------------------------*

!     DETERMINE LOCATIONS OF SUBSTRINGS

 1000 do 200 ifind = 1+offset,linlen
         if (linarr(ifind) .eq. comm(1)) then
            do 205 j = ifind+1,linlen
               if (linarr(j) .eq. comm(2)) then
                  offset = j
                  GOTO 1000
               endif
  205       continue
            GOTO 3000
         endif
         sepchr = .false.
         do 210 isep = 1,nsep
            if ((linarr(ifind) .eq. sep(isep)).or. &
	(linarr(ifind) .eq. tab(isep))) sepchr = .true.
  210    continue
         if (.not. sepchr) then
            nword      = nword + 1
            offset     = ifind
            locb(nword) = ifind
            GOTO 2000
         endif
  200 continue
      GOTO 3000
 2000 do 220 ifind = 1+offset,linlen
         if (linarr(ifind) .eq. comm(1)) then
            loce(nword) = ifind-1
            do 225 j = ifind+1,linlen
               if (linarr(j) .eq. comm(2)) then
                  offset = j
                  GOTO 1000
               endif
225         continue
            GOTO 3000
         endif
         sepchr = .false.
         do 230 isep = 1,nsep
            if ((linarr(ifind) .eq. sep(isep)).or. &
	(linarr(ifind) .eq. tab(isep))) sepchr = .true.
  230    continue
         if (sepchr) then
            loce(nword) = ifind-1
            offset      = ifind
            GOTO 1000
         endif
  220 continue

 3000 continue

!---- EXTRACT SUBSTRINGS -----------------------------------------------

      if (nword .gt. nstrng) GOTO 4100

      do 300 iword = 1,nword
         write(strngs(iword),'(80a)') &
          (linarr(i),i=locb(iword),loce(iword))
  300 continue

!---- JUMP BACK TO CALLING ROUTINE ------------------------------------*

      iret = 0
      errmsg = ' none '
      return

 4100 iret = 1
      errmsg = ' sr PARSE: No. of substrings exceeds NSTRNG* '
      return

      END SUBROUTINE PARSE


! FUNCTION PER LA COMPARAZIONE DI STRINGHE IN MODO CASE-INSENSITIVE
!==================================================================

! ARGOMENTI

! stringa:  la prima stringa (input)
! stringb:  la seconda stringa (input)

! N.B. il riasultato �un logical

	Function conf_strings(stringa,stringb)
	Implicit None
	Logical :: conf_strings
	Character(*), Intent(In) :: stringa,stringb

	Integer :: la,lb,j

	la=len(Trim(AdjustL(stringa)))
	lb=len(Trim(AdjustL(stringb)))
	If(la.Ne.lb) Then
	conf_strings=.False.
	Return
	EndIf

	Do j=1,la
	If((iachar(stringa(j:j)).Ne.iachar(stringb(j:j))).And.&
   	(iachar(stringa(j:j)).Ne.iachar(stringb(j:j))+32).And.&
   	(iachar(stringa(j:j)).Ne.iachar(stringb(j:j))-32)) Then
	conf_strings=.False.
	Return
	EndIf

	conf_strings=.True.

	EndDo

	End Function conf_strings