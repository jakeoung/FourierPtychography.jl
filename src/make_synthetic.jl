using FileIO
using FFTW
using Images

function create_pupil(r, M, N)
    pupil = zeros(M, N)

    for i = floor(Int, M/2-r):ceil(Int, M/2+r)
        for j = floor(Int,N/2-r):ceil(Int,N/2+r)
            if ((i-M/2)^2+(j-N/2)^2<r^2)
                pupil[i,j]=1;
            end
        end
    end
    return pupil
end

function make_synthetic(img_real, phase_real, ratio_LR = 0.1, ratio_step = 0.1*0.4, sigmaN = 0.004)    
    "Return F^{-1}(P * O_l)"
    function get_operator(xx, Masks, n1_LR, n2_LR, pupil)
        L = size(Masks,1)
        xx_c = zeros(ComplexF64, n1_LR, n2_LR, L)
        
        for k = 1:L
            index_x = Masks[k][1]
            index_y = Masks[k][2]

            xx_c[:,:,k] = xx[index_x:(index_x+n1_LR-1),index_y:(index_y+n2_LR-1)] .* pupil
        end
        
        xx_c = ifftshift(xx_c, 1);
        xx_c = ifftshift(xx_c, 2);
        xx_c = ifft( xx_c, [1,2] );
        return xx_c
    end

    img = img_real .* ( cos.(phase_real) .+ 1im*sin.(phase_real) )
    gt = fftshift(fft(img))

    ## 
    H_high, W_high  = size(gt)

    H_low = Int(round(H_high * ratio_LR))
    W_low = Int(round(W_high * ratio_LR))
    step = Int(round(H_high * ratio_step))

    pupil = create_pupil(round(Int, H_low*0.4), H_low, W_low);

    # get mask index
    # k = 0
    kxky = []
    for k1 = round(Int, H_high / 4.0):step:round(Int, H_high/4.0*3)
        for k2 = round(Int, W_high/4):step:round(Int, W_high/4.0*3)
            # k = k + 1;
            push!(kxky, (k1, k2))
        end
    end

    xx_c = get_operator(gt, kxky, H_low, W_low, pupil)
    I = abs.(xx_c) .^ 2

    I_original = copy(I)
    I_original_max = maximum(I_original)
    N_real = sigmaN * maximum(I) * randn(size(I))
    N_real = 0.0
    I = I .+ N_real # I can be negative

    return I, pupil, kxky, img, gt 
end


# img_real = Array{Float64}(load(fimg))
# phase_temp = Array{Float64}(load(fphase))

# phase_temp = (phase_temp .- minimum(phase_temp))
# phase_temp ./= maximum(phase_temp)
# phase_real = phase_temp .* (pi/2.0)

