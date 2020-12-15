program decpot
  !** Calcula la potencia de decaimiento utilizando
  !** los datos de decpow.ini
  use modFlujo
  implicit none

  call leer_ini
  call write_logPower(num,dias,potencias,subcrit,resol)
  if (wrfuen .eq. "si") call write_logSource(num,dias,potencias,subcrit,resol)
  
  call system("PAUSE")

end program decpot
