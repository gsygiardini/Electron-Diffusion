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

function repulsion(board,i,j,p)
    Ly,Lx,Lt = size(board)
        
    if (j < Ly)
        if (board[i,j+1,1] >= 1)
            p[1] = p[1] - 0.1*board[i,j+1,1]
            if (p[1]<0)
                p[1]=0
            end
        end
    end
    if (i < Lx)
        if (board[i+1,j,1] >= 1)
            p[2] = p[2] - 0.1*board[i+1,j,1]
            if (p[2]<0)
                p[2]=0
            end
        end
    end
    if (i > 1)
        if (board[i-1,j,1] >= 1)
            p[3] = p[3] - 0.1*board[i-1,j,1]
            if (p[3]<0)
                p[3]=0
            end
        end
    end
    if (j > 1)
        if (board[i,j-1,1] >= 1)
            p[4] = p[4] - 0.1*board[i,j-1,1]
            if (p[4]<0)
                p[4]=0
            end
        end
    end
    
    sum = p[1]+p[2]+p[3]+p[4]
    p *= (1÷sum) 
    
    return p
end

function process!(board,process)
#     p = [0.25,0.25,0.25,0.25]
    p = [0.4,0.25,0.25,0.1]
    
    Ly,Lx,Lt = size(board)
    R = 0.4*((Lx + Ly)÷4)
    widthI = 0.1  # In fractions of Ly
    widthO = 0.05 # In fractions of Ly
    
    i,j = pop!(process)
    repulsion(board,i,j,p)
    
    move = true
    r = rand()
    ni = i
    nj = j
    #Electron moves forwards
    if (r < p[1])
        if (j < Lx)
            ni = i
            nj = j + 1 
        end
        #Electron moves up
    elseif (r < p[1]+p[2])
        if (i < Ly)
            ni = i + 1
            nj = j
        end
        #Electron moves down
    elseif (r < p[1]+p[2]+p[3])
        if (i > 1)
            ni = i - 1
            nj = j
        end
    else
        #Electron moves backwards
        if (j > 1)
            ni = i
            nj = j - 1
        end
    end
    
        #Electron lifetime dynamics
#       if (board[i,j,2]<=1)
#           board[i,j,1] -= 1
#       end
#       board[ni,nj,2] = board[i,j,2] - 0.1
#       board[i,j,2] = 0
    
    #Channel that leads to the reservoir
    if (nj <= (Lx÷2-R))
        if (ni <= (1-widthI)*(Ly÷2))
            move = false
        elseif (ni >= (1+widthI)*(Ly÷2))
            move = false
        else 
            move = true
        end
    end
    
    #Escape channel
    if (nj >= (Lx÷2+R))
        if (ni <= (1-widthO)*Ly÷2) 
            move = false
        elseif(ni >= (1+widthO)*Ly÷2)
            move = false
        else
            move = true
        end  
    end
    
    #Round"ish" reservoir 
    if (nj > (Lx÷2-R) && nj < (Lx÷2+R))
        if (((ni-Ly÷2)^2 + (nj-Lx÷2)^2) >= R*R)
            move = false
        end
    end
    
    #More than 2 electrons cant occupy the same state
    if (board[ni,nj,1] >= 2)
        move = false
    end

    #if the electron did not move, it does not move
    if(ni==i && nj==j)
        move = false 
    end
    
    #If all conditions are met, the electron can move
    if (move==true)
        #Electrons change after verifying all of the boundary conditions
        board[i,j,1] -= 1
        board[ni,nj,1] += 1
    end
    
    #Delete the electrons that reach the end of the channel
    if (board[i,Lx,1]!=0)
        board[i,Lx,1] = 0
    end
end

function move!(board)
    Ly,Lx,Lt = size(board)
    n = 0
    
    list = []
    for j in 1:Lx
        for i in 1:Ly
            if (board[i,j,1] >= 1)
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
    
    width = 0.1

    j = 1
    for i in (round(Int,(1-width)*Ly÷2)+1):(round(Int,(1+width)*Ly÷2))
        if(board[i,j,1]<2) 
            process = [(i,j)]
            board[i,j,1] += 1
            board[i,j,2] = 10
            process!(board,process) 
        end
    end

    #Altered electrons insertions just for testing
    #i = Ly÷2 + rand(0:floor(Int,Ly÷5)) - Ly÷10
    #j = Lx÷2 + rand(0:floor(Int,Lx÷5)) - Lx÷10

    #i = rand(Ly÷4:(3*Ly)÷4)
    #j = rand(Lx÷4:(3*Lx)÷4)
end

#=========================================================#
#                          MAIN                           #
#=========================================================#
io = open("current.txt", "w")

tMax = 2000
Lx = 50
Ly = 50
Lt = 2
el = 20 # electrons added per timestep
board = zeros(Ly,Lx,Lt)
board[:,:,2] .= 0.0

#Vector for counting the number of electrons that go through
current = zeros(Int64,tMax)

print("Starting simulation...\n")
# if (false)
anim = @animate for k = 1:tMax
    heatmap(board[:,:,1])
    action!(board)
    move!(board)
    
    if (mod(k,100)==0)
        print("$k of $tMax \n")
    end
    print(io,"$k " * string(count(board)) * " \n")
    current[k]=count(board)    
end

#Current Fourier transform - not ready yet
#t = 1:tMax
#R = real(current)
# R2 *= R
#I = imag(current)
# I2 *= I

#plot(1/t,R,title="Current")
#savefig("myplot.png")

gif(anim, "test.gif", fps = 20)
# end
