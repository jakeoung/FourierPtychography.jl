using Plots

"Gerchberg-Saxton method"
function recon_gs(I, O0, kxky, pupil, niter, dresult, gt=nothing)
    mkpath(dresult)
    O = copy(O0)

    H_high, W_high = size(O0, 1), size(O0, 2)
    H_low,  W_low   = size(I, 1), size(I, 2)

    pupil_bool = (pupil .>= 0.1)

    for loop = 1 : niter
        for num = 1 : length(kxky)
            # Get the subspectrum
            index_x = round(Int, kxky[num][1])
            index_y = round(Int, kxky[num][2])

            Subspectrum = O[index_x:(index_x+H_low-1), index_y:(index_y+W_low-1)]
            
            # Subspectrum = F[Fcenter_Y+ky-Int(trunc(M/2)): Fcenter_Y+ky+Int(ceil(M/2))-1, Fcenter_X+kx-Int(trunc(N/2)):Fcenter_X+kx+Int(ceil(M/2))-1]
            
            Abbr_Subspectrum = Subspectrum .* pupil # P(u) .* O(u)
            
            # Real space modulus constraint
            # Uold = ifft(fftshift(Abbr_Subspectrum))
            # Abbr_Subspectrum = ifftshift(Abbr_Subspectrum, 1)
            # Abbr_Subspectrum = ifftshift(Abbr_Subspectrum, 2)
            Abbr_Subspectrum = ifftshift(Abbr_Subspectrum)
            Uold = ifft(Abbr_Subspectrum)
            Unew = I[:,:,num] .* (Uold./abs.(Uold))
            
            # Fourier space constraint and object function update
            Abbr_Subspectrum_corrected = fftshift(fft(Unew))
            
            W = abs.(pupil)./maximum(abs.(pupil))
            
            invP = conj(pupil)./((abs.(pupil)) .^2 .+ eps() .^ 2)
            Subspectrumnew = (W.*Abbr_Subspectrum_corrected .+ (1 .- W) .* (Abbr_Subspectrum)) .* invP
            Subspectrumnew[pupil_bool.==0] = Subspectrum[pupil_bool.==0]
            
            # Fourier spertrum replacement
            O[index_x:(index_x+H_low-1), index_y:(index_y+W_low-1)] = Subspectrumnew
            # F[Fcenter_Y+ky-Int(trunc(M/2)): Fcenter_Y+ky+Int(ceil(M/2))-1, Fcenter_X+kx-Int(trunc((N/2))):Fcenter_X+kx+Int(ceil(M/2))-1] = Subspectrumnew
            
            # Inverse FT to get the reconstruction
            # if num == LED_num_x*LED_num_y
            #     Result .= ifft(fftshift(F))
            # end
        end
        if loop % 1 == 0
            im_abs = abs.(ifft(fftshift(O)))
            # p1 = heatmap(log.(abs.(O).+1), title="Fourier spectrum")
            # p2 = heatmap(abs.(ifft(fftshift(O))), title="Reconstructed amplitude $num")
            # save(dresult*"$t.png", im_r ./ maximum(im_r))
            # plot(p1, p2)
            save(dresult*"$(loop).png", im_abs ./ maximum(im_abs))
        end
    end
    return O
end