MODULE Task
use fgsl
use, intrinsic :: iso_c_binding
implicit none
    integer( fgsl_size_t ), private    :: n      !Число узлов.
    type( fgsl_interp_accel ), private :: acc    !Акселератор.
    type( fgsl_spline ), private       :: spline !Кубический сплайн.
    integer( fgsl_int ), private       :: st     !Статус.
contains


    !Выполняет интерполяцию кубическим сплайном.
    FUNCTION GetSplineValue( x, x_array, y_array )
    real( fgsl_double ), dimension(:), intent(in) :: x_array        !Узлы интерполяции.
    real( fgsl_double ), dimension(:), intent(in) :: y_array        !Значения функции в узлах.
    real( fgsl_double ), intent(in)               :: x              !Точка, в которой ищем значение.
    real( fgsl_double )                           :: GetSplineValue !Значение f(x).

        n = size( x_array, dim = 1 )                          !Получает число узлов.
        spline = fgsl_spline_alloc( fgsl_interp_cspline, n )  !Выделяет память под сплайн.
        acc = fgsl_interp_accel_alloc()                       !Выделяет память под акселератор.
        st = fgsl_spline_init( spline, x_array, y_array )     !Инициализация сплайна.
        GetSplineValue = fgsl_spline_eval( spline, x, acc )   !Получает f(x).Определяет сплайн.

        !Очистка памяти.
        call fgsl_spline_free( spline )
        call fgsl_interp_accel_free( acc )

    END FUNCTION


    !Задаем функцию для интегрирования.
    FUNCTION f( x, params ) bind( c )
    real( c_double ), value   :: x
    real( c_double )          :: f
    type( c_ptr    ), value   :: params
    real( c_double ), pointer :: arg

        call c_f_pointer( params, arg )
        f = fgsl_spline_eval( spline, x, acc )

    END FUNCTION


    !Выполнете интегрирование сплайна.
    FUNCTION GetSplineIntegral( a, b, x_array, y_array )
    real( fgsl_double ), dimension(:), intent(in) :: x_array           !Узлы интерполяции.
    real( fgsl_double ), dimension(:), intent(in) :: y_array           !Значения функции в узлах.
    real( fgsl_double ), intent(in)               :: a                 !Верхний предел интегрирования.
    real( fgsl_double ), intent(in)               :: b                 !Нижний предел интегрирования.
    real( fgsl_double )                           :: GetSplineIntegral !Значение интеграла.
    real( fgsl_double), parameter                 :: epsabs = 0.0
    real( fgsl_double), parameter                 :: epsrel = 1.0e-7
    integer( fgsl_size_t ), parameter             :: nmax = 1000
    real( fgsl_double )                           :: errest            !Оценка ошибки интегрирования.
    real( fgsl_double ), target                   :: arg
    type( c_ptr )                                 :: ptr
    type( fgsl_function )                         :: f_obj
    type( fgsl_integration_workspace )            :: wk

        !Получение сплайна.
        n = size( x_array, dim = 1 )                          !Получает число узлов.
        spline = fgsl_spline_alloc( fgsl_interp_cspline, n )  !Выделяет память под сплайн.
        acc = fgsl_interp_accel_alloc()                       !Выделяет память под акселератор.
        st = fgsl_spline_init( spline, x_array, y_array )     !Инициализация сплайна.

        ptr = c_loc( arg )

        !Интегрирование сплайна.
        f_obj = fgsl_function_init( f, ptr )
        wk = fgsl_integration_workspace_alloc( nmax )
        st = fgsl_integration_qags( f_obj, a, b, epsabs, epsrel, nmax, wk, GetSplineIntegral, errest ) !Само интегрирование.

        !Очистка памяти.
        call fgsl_integration_workspace_free( wk )
        call fgsl_function_free( f_obj )
        call fgsl_spline_free( spline )
        call fgsl_interp_accel_free( acc )

    END FUNCTION


END MODULE
