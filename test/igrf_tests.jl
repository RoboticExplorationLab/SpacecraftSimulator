using Test, Statistics


@testset "igrf 13" begin
    let

        date = 2020.3
        lat = 34.567
        elong = 45.678
        alt = 1.05*R_EARTH/1000

        B_true1 = [24354.4;1908.4;32051.6]
        B1   = my_igrf_13(date,alt,lat,elong,13)

        @test isapprox(B1,B_true1,atol=1e-1)

        # test case two
        date = 2024.6
        lat = -78.123
        elong = 150.4
        alt = 1.10*R_EARTH/1000

        B_true2 = [-6554.1;1913.3;-44845.9]
        B2  = my_igrf_13(date,alt,lat,elong,13)

        @test isapprox(B2,B_true2,atol=1e-1)

        # test case three
        date = 2024.5
        lat = -35.78
        elong = -130.48
        alt = 1.07*R_EARTH/1000

        B_true3 = [19285.7;6911.3;-25084.9]
        B3  = my_igrf_13(date,alt,lat,elong,13)

        @test isapprox(B3,B_true3,atol=1e-1)

    end
end

@testset "5th order IGRF" begin
    let


        # test the 5th order one
        date = 2020.3
        lat = 34.567
        elong = 45.678
        alt = 1.05*R_EARTH/1000

        B13_5_1   = my_igrf_13(date,alt,lat,elong,5)
        B5_5_1   = my_igrf_5(date,alt,lat,elong)

        @test isapprox(norm(B13_5_1 - B5_5_1),0,rtol = 1e-9)

        # test 5th order one again
        date = 2020.3
        lat = 34.567
        elong = 45.678
        alt = 1.05*R_EARTH/1000

        B13_5_1   = my_igrf_13(date,alt,lat,elong,5)
        B5_5_1   = my_igrf_5(date,alt,lat,elong)

        @test isapprox(norm(B13_5_1 - B5_5_1),0,rtol = 1e-9)

    end
end

mag_err = zeros(1000)
ang_err = zeros(1000)
@testset "errors for 5th order IGRF" begin
    let

        for i = 1:1000
            date = rand_in_range(2020,2024)
            lat = rand_in_range(-90,90)
            elong = rand_in_range(-180,180)
            alt = rand_in_range(1,1.1)*R_EARTH

            B_13  = my_igrf_13(date,alt,lat,elong,13)
            B_5   = my_igrf_5(date,alt,lat,elong)

            # check angle and attitude errors
            ang_err[i] = acos(dot(normalize(B_13),normalize(B_5)))
            mag_err[i] = 100*(norm(B_5)-norm(B_13))/norm(B_13)

        end

        # test max
        @test maximum(rad2deg.(ang_err))<10.0
        @test mean(rad2deg.(ang_err))<3.0

        @test maximum(mag_err) < 10
        @test mean(mag_err) < 5

    end
end
