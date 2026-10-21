using Test
include("../utils/generar_csv.jl")
include("../utils/resultador.jl")

ENV["GKSwstype"] = "100"

@testset "Comparación genérica de dos métodos" begin
    mktempdir() do dir
        t, S = [0.0, 0.1], ComplexF64[1, im]
        llamadas_a, llamadas_b, llamadas_ref = Ref(0), Ref(0), Ref(0)
        a = N -> begin
            llamadas_a[] += 1
            t, S .* (1 + (1,3,5)[mod1(llamadas_a[],3)]/(100N))
        end
        b = N -> begin
            llamadas_b[] += 1
            t, S .* (1 + 1/(100N))
        end
        ref = () -> (llamadas_ref[] += 1; (t,S))
        path = joinpath(dir,"errores.csv")
        df = generar_csv_errores(path,a,b,ref; nombres=("Aleatorio","Determinista"),
            repeticiones=(3,1), Ns=[6,2,4,1], verbose=false)
        @test (llamadas_a[], llamadas_b[], llamadas_ref[]) == (12,4,1)
        @test nrow(df) == 8
        A, B = df[df.method.=="Aleatorio",:], df[df.method.=="Determinista",:]
        @test A.error_mean_pct ≈ 3 ./ A.N_iso
        @test A.error_sd_pct ≈ 2 ./ A.N_iso
        @test A.error_sem_pct ≈ 2 ./ (A.N_iso .* sqrt(3))
        @test all(ismissing,B.error_sem_pct) && all(ismissing,B.error_sd_pct)
        @test isequal(leer_errores_csv(path), df)
        @test nrmse_complejo_porcentual(S,im .*S) ≈ 100sqrt(2)
        cruce = buscar_equivalente_csv(df,6;verbose=false)
        @test cruce.N_equivalente == 2
        @test cruce.N_interpolado ≈ 2
        @test cruce.reduccion_pct ≈ 100*(1-2/6)
        @test cruce.factor == 3
        inverso = buscar_equivalente_csv(df,1;metodos=("Determinista","Aleatorio"),verbose=false)
        @test inverso.N_equivalente == 4
        p, resumen = graficar_errores_csv(path;objetivos=[6],mostrar=false,verbose=false)
        @test resumen.N_equivalente == [2]
        @test p[1][:yaxis][:scale] == :log10
        @test filesize(joinpath(dir,"errores_equivalencias.png")) > 0
        @test filesize(joinpath(dir,"errores_equivalencias.pdf")) > 0
        @test p.series_list[1][:ribbon][2] ≈ 1.96 .* A.error_sem_pct

        # Ceros: escala lineal. Dos deterministas: sin bandas y sin equivalencias obligatorias.
        generar_csv_errores(path,N->(t,S),N->(t,-S),(t,S);verbose=false)
        p, resumen = graficar_errores_csv(path;guardar=false,mostrar=false,verbose=false)
        @test p[1][:yaxis][:scale] == :identity
        @test length(p.series_list) == 2
        @test nrow(resumen) == 0

        # Ambos pueden tener repeticiones; el CSV se reproduce al fijar la semilla.
        aleatorio = N -> (t,S .* (1+rand()/N))
        d1 = generar_csv_errores(path,aleatorio,aleatorio,(t,S);
            repeticiones=(4,3),Ns=2:4,seed=19,verbose=false)
        d2 = generar_csv_errores(path,aleatorio,aleatorio,(t,S);
            repeticiones=(4,3),Ns=2:4,seed=19,verbose=false)
        @test isequal(d1, d2)
        @test all(!ismissing, d1.error_sem_pct)
        @test_throws ArgumentError generar_csv_errores(path,N->(t.+.01,S),b,(t,S);verbose=false)
        @test_throws ArgumentError generar_csv_errores(path,b,b,(t,S);repeticiones=(0,1),verbose=false)
        @test_throws ArgumentError generar_csv_errores(path,b,b,(t,S);nombres=("X","X"),verbose=false)
        bad = copy(d1); bad.error_metric .= "complex_mse"
        @test_throws ArgumentError validar_errores_csv(bad)
        bad = vcat(d1,d1[1:1,:])
        @test_throws ArgumentError validar_errores_csv(bad)

        # Sin cruce y cruce con error cero: no requieren log(0).
        d = generar_csv_errores(path,N->(t,S),b,(t,S);Ns=1:4,verbose=false)
        @test ismissing(buscar_equivalente_csv(d,4;verbose=false).N_equivalente)
        d = generar_csv_errores(path,b,N->(t,N<3 ? 2S : S),(t,S);Ns=1:4,verbose=false)
        @test buscar_equivalente_csv(d,4;verbose=false).N_equivalente == 3
    end
end

@testset "CSV porcentual anterior" begin
    old = joinpath(@__DIR__,"..","errores_FID.csv")
    if isfile(old)
        original = CSV.read(old,DataFrame)
        adapted = leer_errores_csv(old)
        @test nrow(adapted) == 2nrow(original)
        @test adapted.error_mean_pct[adapted.method.=="Monte Carlo"] == original.mc_mean_pct
        @test adapted.error_mean_pct[adapted.method.=="QMC-Halton"] == original.halton_error_pct
    end
end
