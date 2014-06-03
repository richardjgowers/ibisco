!> @file
!> @brief Helping variables for parsing

MODULE moduleparsing

 Character(80) :: line,errmesg
 Integer, Parameter :: nsep=1,nstrng=50
 Character(80) :: strngs(nstrng)
 Character(1) :: sep(nsep),comm(2), tab(nsep)
 Integer :: nword,iret

data sep /' '/
data tab /'	'/
data comm(1) /'#'/
data comm(2) /'#'/

END MODULE moduleparsing
