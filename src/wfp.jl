function get_operator(xx, kxky, n1_LR, n2_LR, pupil)
    
    L = size(kxky,1)
    xx_c = zeros(ComplexF64, n1_LR, n2_LR, L)
    
    for k = 1:L
        index_x = kxky[k][1]
        index_y = kxky[k][2]

        xx_c[:,:,k] = xx[index_x:(index_x+n1_LR-1),index_y:(index_y+n2_LR-1)] .* pupil
    end
    
    xx_c = ifftshift(xx_c, 1);
    xx_c = ifftshift(xx_c, 2);
    xx_c = ifft( xx_c, [1,2] );
    return xx_c
end

function get_inv_operator(xx_c, kxky, n1, n2, pupil)
    n1_LR, n2_LR, L = size(xx_c);

    xx_c = fft(xx_c, (1,2))
    xx_c = fftshift(xx_c,1)
    xx_c = fftshift(xx_c,2)

    min1 = minimum(map(i -> i[1], kxky))
    max1 = maximum(map(i -> i[1], kxky))
    min2 = minimum(map(i -> i[2], kxky))
    max2 = maximum(map(i -> i[2], kxky))
    
    xx = zeros(ComplexF64, max1 - min1  + n1_LR, max2 - min2 + n2_LR, L)

    for k = 1:L
        index_x = kxky[k][1];
        index_y = kxky[k][2];
        xx[ index_x-min1+1 : index_x-min1+n1_LR , index_y-min2+1 : index_y-min2+n2_LR , k ] .= xx_c[:,:,k].*conj(pupil)
    end

    xx_mean = dropdims(sum(xx, dims=3), dims=3) / size(xx, 3)

    xx_mean_large = zeros(ComplexF64, n1,n2)
    xx_mean_large[ min1:min1+size(xx_mean,1)-1 , min2:min2+size(xx_mean,2)-1 ] .= xx_mean

    return xx_mean_large
end

"""
    recon_wfp

# Args
- I [H_high x W_high x L] : measurements
- O0 : initial guess on the Fourier domain
- kxky : Mask indices

# Returns
- O : 

The code is ported to Julia from https://www.sites.google.com/site/lihengbian/
"""
function recon_wfp(I, O0, kxky, pupil, niter, mu_max, sigma2, weight, dresult, gt=nothing)
    H_high, W_high = size(O0, 1), size(O0, 2)
    H_high_LR, W_high_LR, L = size(I)

    O = copy(O0)

    Relerrs = zeros(niter+1)
    (~isnothing(gt)) && (Relerrs[1] = sum( abs.(abs.(O) - abs.(gt)) ) / sum(abs.(gt)))
    # Relerrs( 1 ) = sum(sum(abs(abs(z)-abs(gt))))/sum(sum(abs(gt))) ;
    N = zeros(size(I));
    epsilon = N;

    normest = sqrt(sum(I)/ L / H_high_LR / W_high_LR)
    StepsizeN = 0.01;

    im_r = zeros(size(O))

    #O is in the F domain
    ## iter
    for t=1:niter
        @show t
        # update O
        Amplitude_r = get_operator(O, kxky, H_high_LR, W_high_LR, pupil)
        gf = ((abs.(Amplitude_r)).^2 + N - I )  .* Amplitude_r
        w = get_inv_operator(gf, kxky, H_high, W_high, pupil)
        
        alpha = - log(0.997)
        mu  = 1-exp(-alpha*t)
        mu  = min(mu, mu_max) / (normest^2)
        # mu *= 0.01
        O   = O - mu * w
        
        # update N
        CN = (abs.(Amplitude_r).^2 + N - I ) + weight .* (N.*N .- 9*sigma2 .+ epsilon.*epsilon) .*(2*N)
        stepN = mu * StepsizeN
        N = N - stepN * CN
        
        # update epsilon
        Etemp = 9*sigma2 .- N.*N
        Etemp[Etemp.<0] .= 0
        epsilon = sqrt.(Etemp)
        
        #######################################################################
        # calculate recovery error
        (~isnothing(gt)) && (Relerrs[t + 1] = sum(abs.(abs.(O)-abs.(gt))) / sum(abs.(gt)))
        # save iterative results
        if mod(t,50) == 0
            im_r = ifft(ifftshift(O))
            abs_im_r = abs.(im_r)
            save(dresult*"$t.png", abs_im_r ./ maximum(abs_im_r))
            im_r_ang = angle.(im_r)/pi;
            im_r_ang = im_r_ang .- minimum(im_r_ang)
            im_r_ang = abs.(im_r_ang ./ maximum(im_r_ang))
            save(dresult*"$(t)_ang.png", im_r_ang ./ maximum(im_r_ang))  
        end
        #######################################################################
    end
   
    return O, im_r, Relerrs
end
