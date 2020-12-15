module modFlujo
  !******************************************************************!
  ! Sencillo módulo con utilidades para el cálculo de la potencia 
  ! de decaimiento de un reactor teniendo en cuenta la historia de 
  ! potencia, las fracciones de precursores y sus lambdas
  implicit none
  integer :: num_grupos = 0     !Cantidad de grupos de precursores.
  integer :: num = 0            !numero de períodos (potencias)
  integer :: resol              !Número de puntos totales a calcular
  character(len=13) :: archivo = "decpow.ini" !Archivo inicial con datos
  character(len=2) :: wrfuen = "no"      !Escribir los resultados para la fuente también?        !
  real(kind(0d0)) :: subcrit             !Subcriticidad del núcleo
  real(kind(0d0)), allocatable :: beta(:), lambda(:) 
  real(kind(0d0)), allocatable :: dias(:), potencias(:) 
  contains

  subroutine leer_ini
  !* Intenta leer el archivo del configuración
    implicit none
    integer :: i
    character(len=5) dumm, rfuen*1 
    
    open(55, file=archivo, status='old')
    write(*,*)
    write(*,'(A40)') "############# DECPOW ################"
    write(*,'(A8,A13)') "Leyendo ",archivo

    read(55,*) !Resolución
    read(55,'(I10)') resol
    write(*,'(A12,I10,A7)') "Resolucion: ", resol, " puntos"  

    read(55,*) !Subcriticidad
    read(55,'(F9.3)') subcrit
    write(*,'(A15,F7.3,A3)') "Subcriticidad: ", subcrit, ' mk'

    read(55,*) !Escribir fuente?
    read(55,'(A1)') rfuen
    if ((rfuen .eq. "S") .or. (rfuen .eq. "s")) wrfuen="si"
    write(*,'(A32,A2)') "Escribir archivo de la fuente?: ", wrfuen

    read(55,*) !grupos de fotoneutrones
    i = 0
    leer_foton: do
      i=i+1
      read(55,'(A5)') dumm
      if (dumm .eq. "TRANS") exit leer_foton
    enddo leer_foton
    num_grupos=i-1
    write(*,'(A25,I2)') "Grupos de fotoneutrones: ", num_grupos
 
    !Transitorio
    i = 0
    leer_transitorio: do
      i=i+1
      read(55,*) dumm
      if ( trim(adjustl(dumm)) .eq. "FIN") exit leer_transitorio
    enddo leer_transitorio
    num=i-1
    write(*,'(A15,I3,A11)') "Transitorio de ", num, " potencias."
    close(55)

    allocate(beta(num_grupos)); allocate(lambda(num_grupos))
    allocate(dias(num)); allocate(potencias(num))
    open(55,file="candecpow.ini", status='old')
    read(55,'(6/)')
    do i=1,num_grupos
      read(55,*)beta(i), lambda(i)
    enddo
    read(55,*) !TRANS
    do i=1,num
      read(55,*)dias(i),potencias(i)
    enddo
    close(55)
    write(*,'(A19,F9.3,A6)') "Intervalo total de ", sum(dias), " dias."
    write(*,'(A40)') "############# DECPOW ################"
    write(*,*) 
  end subroutine leer_ini

  subroutine get_logPower(N,days,powers,subcrit_mk,resolution,time,logPower)
    !! Calcula la potencia total, como fracción de PP
    !! Guarda los resultados en memoria, permite seguir con otros calculos, 
    !!pero puede consumir mucha RAM si la resolucion es demasiado alta.
    implicit none
    integer, intent(in) :: N, resolution
    real(kind(0d0)), intent(in) :: days(N), powers(N), subcrit_mk
    real(kind(0d0)), intent(out) :: time(resolution), logPower(resolution)
    real(kind(0d0)) :: periods(N+1)
    integer :: i,j

    call create_periods(N,days,periods)

    do i = 1, resolution
      time(i) = periods(N+1)*real(i,kind(0d0))/real(resolution,kind(0d0))
      do j=1,N
        if ( (time(i) .ge. periods(j)) .and. (time(i) .le. periods(j+1)) ) then
          logPower(i) =  max(dlog10(sourceFFP(subcrit_mk,total_source(N,periods,powers,time(i)))),dlog10(powers(j)))
        endif
      enddo
    enddo
  end subroutine get_logPower

  subroutine write_logPower(N,days,powers,subcrit_mk,resolution)
    !! Calcula la potencia total, como fracción de PP.
    !! Escribe directamente los resultados, sin guardarlos en memoria, 
    !! por si la RAM no alcanza.
    implicit none
    integer, intent(in) :: N, resolution
    real(kind(0d0)), intent(in) :: days(N), powers(N), subcrit_mk
    real(kind(0d0)) :: periods(N+1), time
    integer :: i,j

    call create_periods(N,days,periods)

    open(59,file='logPower.dat')
    do i = 1, resolution
      time = periods(N+1)*real(i,kind(0d0))/real(resolution,kind(0d0))
      do j=1,N
        if ( (time .ge. periods(j)) .and. (time .le. periods(j+1)) ) then
          write(59,*)time/86400.d0, max(dlog10(sourceFFP(subcrit_mk,total_source(N,periods,powers,time))),dlog10(powers(j)))
        endif
      enddo
    enddo
    close(59)
  end subroutine write_logPower

  subroutine get_logSource(N,days,powers,subcrit_mk,resolution,time,logSource)
    !! Calcula la potencia de la fuente, como fracción de PP
    !! Guarda los resultados en memoria, permite seguir con otros calculos, 
    !!pero puede consumir mucha RAM si la resolucion es demasiado alta.
    implicit none
    integer, intent(in) :: N, resolution
    real(kind(0d0)), intent(in) :: days(N), powers(N), subcrit_mk
    real(kind(0d0)), intent(out) :: time(resolution), logSource(resolution)
    real(kind(0d0)) :: periods(N+1)
    integer :: i

    call create_periods(N,days,periods)

    do i = 1, resolution
      time(i) = periods(N+1)*real(i,kind(0d0))/real(resolution,kind(0d0))
      logSource(i) =  dlog10(sourceFFP(subcrit_mk,total_source(N,periods,powers,time(i))))
    enddo
  end subroutine get_logSource

  subroutine write_logSource(N,days,powers,subcrit_mk,resolution)  !!VER hacer optional la resolucion. Falta definir View!!
    !! Calcula la potencia de la fuente, como fracción de PP.
    !! Escribe directamente los resultados, sin guardarlos en memoria, 
    !! por si la RAM no alcanza.
    implicit none
    integer, intent(in) :: N, resolution
    real(kind(0d0)), intent(in) :: days(N), powers(N), subcrit_mk
    real(kind(0d0)) :: periods(N+1), time
    integer :: i

    call create_periods(N,days,periods)

    open(58,file='logSource.dat')
    do i = 1, resolution
      time = periods(N+1)*real(i,kind(0d0))/real(resolution,kind(0d0))
      write(58,*)time/86400.d0, dlog10(sourceFFP(subcrit_mk,total_source(N,periods,powers,time)))
    enddo
    close(58)
  end subroutine write_logSource

  subroutine create_periods(N,days,periods)
    implicit none
    integer, intent(in) :: N
    real(kind(0d0)), intent(in) :: days(N)             !! Dias
    real(kind(0d0)), intent(out) :: periods(N+1)       !! tiempos límite, en segundos
    integer :: i
    periods(1)=0.d0
    do i=2,N+1
      periods(i) = periods(i-1)+days(i-1)*24.d0*60.d0*60.d0
    enddo
  end subroutine create_periods

  real(kind(0d0)) function sourceFFP(subcrit_mk,source)
    implicit none
    real(kind(0d0)) :: subcrit_mk, source
    sourceFFP = 1000.d0*source/subcrit_mk
  end function sourceFFP

  real(kind(0d0)) function total_source(N, periods, powers, time)
    implicit none
    integer :: N
    real(kind(0d0)) :: periods(N+1), powers(N), time                                !! Periods posee los tiempos de corte 0,t1,t2,...tN
    real(kind(0d0)) :: build, decay                                                 !! powers posee las potencias P1,P2,...PN
    integer :: i
    total_source = 0.d0
    do_loop: do i = 1, N
      build = min(periods(i+1),time) - periods(i)                                   ![t_i,min(t_i+1,time)]
      decay = max(0.d0, time - periods(i+1))                                        ![t_i+1,time]
      total_source = total_source + partial_source(build,decay,powers(i))
      if (time < periods(i+1)) exit do_loop                                         !Si t < t_i+1, el resto de los terminos es nulo
      enddo do_loop 
  end function total_source

  real(kind(0d0)) function partial_source(buildTime,decayTime,power)
    implicit none
    !Construye la fuente, para una potencia, tiempo de buildup y un tiempo de decaimiento dados (en segundos).
    real(kind(0d0)) :: buildTime, decayTime, power
    integer :: i
    partial_source = 0.d0
    do i=1, num_grupos
      partial_source = partial_source + 1.d-2 * beta(i) * power * buildUp(lambda(i),buildTime) * decay(lambda(i),decayTime)
    enddo
  end function partial_source

  real(kind(0d0)) function buildUp(lambda,time)
    !* Factor de buildUp, para un lambda y un tiempo dado
    implicit none
    real(kind(0d0)) :: lambda, time
    buildUp = 1 - dexp(-lambda*time)
  end function buildUp

  real(kind(0d0)) function decay(lambda,time)
    !* Factor de decaimiento exponencial, para un lambda y un tiempo dado
    implicit none
    real(kind(0d0)) :: lambda, time
    decay = dexp(-lambda*time)
  end function decay

end module modFlujo
