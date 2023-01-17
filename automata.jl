using Plots
using FFTW

function count(board)
    Ly,Lx,Lt = size(board)
    
    n=0
    j=Lx-1
    for i in 1:Ly
        if (board[i,j,1]>0)
            n += 1 
        end
    end
    return n
end

function critical(board,y,x)
    Ly,Lx,Lt = size(board)
    c = false

    if (y > 1)
        if (board[y-1,x,1] == 1)
            c = c | true
        end
    end
    if (y < Ly)
        if (board[y+1,x,1] == 1)
            c = c | true
        end
    end
    if (x > 1)
        if (board[y,x-1,1] == 1)
            c = c | true
        end
    end
    if (x < Lx)
        if (board[y,x+1,1] == 1)
            c = c | true
        end
    end

    return c
end

function process!(board,process)
    p1 = 0.40
    p2 = 0.60
    p3 = 0.80
    
    Ly,Lx,Lt = size(board)
    R = 0.5*(Lx + Ly)÷4
    
    moved = false

    while(length(process) > 0)
#     while (moved == false)
        i,j = pop!(process)
        place = false
        r = rand()
        ni = i
        nj = j
        if (r < p1)
            if (j < Lx)
                ni = i
                nj = j + 1
                place = true
            end
        elseif (r < p2)
            if (j > 1)
                ni = i
                nj = j - 1
                place = true
            end
        elseif (r < p3)
            if (i < Ly)
                ni = i + 1
                nj = j
                place = true
            end
        else
            if (i > 1)
                ni = i - 1
                nj = j
                place = true
            end
        end
        
        #Electron deleted before jumping
#         board[i,j,1] -= 1
        
        #Electron lifetime dynamics
#       if (board[i,j,2]<=1)
#           place = false
#       end
#       board[ni,nj,2] = board[i,j,2] - 0.1
    #   board[ni,nj,2] = board[i,j,2]
    #   board[i,j,2] = 0
            
        if (place == true)
            #Electron changing neighboors before boundary conditions
#             board[ni,nj,1] += 1
#             board[i,j,1] -= 1

            c = board[ni,nj,1] > 1   
#             c |= board[ni,nj,1] ≥ 1   
            #Fermionic
            c |= critical(board,ni,nj)                                          # e-e repulsion
            c |= (((nj > 2 && nj <= 4) || (nj > 6 && nj <= 8)) && ni > 6)
#            dot
            c |= ni ≤ 1 || ni ≥ Ly                                             # walls
            c |= nj ≤ 1 || nj ≥ Lx                                             # walls
            c |= ((ni-Ly÷2)^2 + (nj-Lx÷2)^2) ≥ R*R
#             round walls
            
            if (c)
                #Electrons change after verifying all of the doundary conditions
                board[i,j,1] -= 1
                board[ni,nj,1] += 1
                push!(process,(ni,nj))
                
                #Condition that forces electron to change place (may stuck the program if there is no available neighboor to move)
#                 moved=true
            end
        end
    end
#     end
end

function move!(board)
    Ly,Lx,Lt = size(board)
    n = 0
    
    list = []
    for j in 1:Lx
        for i in 1:Ly
            if (board[i,j,1] == 1)
                push!(list,(i,j))
                n+=1
            end
        end
    end
    #Instead of randomly choosing a position and moving the electron
    #1. A list is created with all non zero positions 
    #2. A random element from the list is selected and moved
    #3. The moved electron is excluded from the list such that
    #each electron is moved only once per time step
    #Otherwise, we observe electrons that do not move and electrons that move multiple times
    for i in 1:n
        r = rand(1:length(list))
        i,j = list[r]
        list = deleteat!(list,r)
    
        process = [(i,j)]
        process!(board,process)
    end
end

function action!(board)
    Ly,Lx,Lt = size(board)
#     i = Ly÷2
#     j = 1
#     j = Lx÷2
    #Altered electrons insertions just for testing
    i = rand(Ly÷4:(3*Ly)÷4)
    j = rand(Lx÷4:(3*Lx)÷4)

    process = [(i,j)]
    board[i,j,1] += 1
    board[i,j,2] = 10

    process!(board,process)
end

#=========================================================#
#                          MAIN                           #
#=========================================================#
io = open("current.txt", "w")

tMax = 20
Lx = 50
Ly = 50
Lt = 2
el = 100
board = zeros(Ly,Lx,Lt)
board[:,:,2] .= 0.0

# current = zeros(Int64,tMax,2)
current = zeros(Int64,tMax)

#Function for adding electrons at the beggining of the simulation just for testing
for i in 1:el
    action!(board)
end

# if (false)
anim = @animate for k = 1:tMax
    heatmap(board[:,:,1])
#     action!(board)

    move!(board)
    
    print("$k \n")
    print(io,"$k " * string(count(board)) * " \n")
    current[k]=count(board)    
end

#Current Fourier transform - not ready yet
# t = 1:1/tMax:tMax
# plot(t,real(current),title="Current")

# plot(t,(real(fft(current))*real(fft(current)) + imag(fft(current))*imag(fft(current))),title="Signal")
# savefig("myplot.png")

gif(anim, "test.gif", fps = 1)
# end
