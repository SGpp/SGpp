from matplotlib import pyplot as plt
import numpy as np
import pysgpp


def testplot():

    X = np.linspace(0, 1, 900)

    psiBase = pysgpp.SPsiHermiteBase()
    zetaBase = pysgpp.SZetaHermiteBase()

    psi_values = [psiBase.eval(1, 1, x) for x in X]
    zeta_values = [zetaBase.eval(1, 1, x) for x in X]

    return X, psi_values, zeta_values



def calc_y_values(X,function,level):
    X_all=[]
    Y_all=[]
    for i in range( len(X)):



        X_samples,Y_samples=getrelevantValues(X[i][0],function,X[i][1],X[i][2])
        X_all.append(X_samples)
        Y_all.append(Y_samples)


    return X_all,Y_all





def getrelevantValues(x,function,level,index):
    h_l=2**-level
    X = np.linspace(x-h_l, x+h_l, 300)




    Y=[trans_dilatate(level,index,point,function) for point in X]

    return X,Y





def psi(x):
    if (x < 0):
        return -2 * x ** 3 - 3 * x ** 2 + 1
    else:
        return 2 * x ** 3 - 3 * x ** 2 + 1


def zeta(x):
    if(x < 0):
        return x ** 3 + 2 * x ** 2 + x
    if(x >= 0):
        return x ** 3 - 2 * x ** 2 + x


def hat(x):
    return  max(1-abs(x),0)



def createGrid(level):
    if(level==0):

        return [(0.5,level,1)]
    #elif(level==1):
    #    return [(0,level,0),(1,level,2)]

    else:
        h_l = 2 ** -(level)
        grid=[]

        for i in range((2**level)+1):
            #if(i%2==0):
            #    continue

            triple = (i * h_l,level,i)
            grid.append(triple)


        return grid



def trans_dilatate(level, index,x,function):
    if(level==0):
        h_l=0.5
        index=1

    else:


        h_l=2**-level





    return function(x/h_l-index)


def uniformZetaPsi():
    X = np.linspace(-1, 1, 900)
    psi_values = [psi(x) for x in X]
    zeta_values = [zeta(x) for x in X]

    return X, psi_values, zeta_values


#X, psi_values, zeta_values = testplot()

level=4


def linearnodalgrid():
    f, axarr = plt.subplots(4, 1)

    ax = axarr[1]

    for l in range(level):
        X = createGrid(l)

        X_grids, Y_values = calc_y_values(X, hat, l)
        ax = axarr[l]

        for i in range(len(X_grids)):
            if (X[i][2] == 0 or X[i][2] == len(X_grids) - 1):
            #if(X[i][2]%2==0 and X[i][1]!=0):
               ax.plot(X_grids[i], Y_values[i], linestyle="--", c="grey")
            else:
                ax.plot(X_grids[i], Y_values[i])

        # ax.plot(X,  psi_values, label='$\psi$')
        # ax.plot(X, zeta_values, label="$\zeta$")
        # plt.legend()

        # x axis "labels"
        ax.set_xticks([x[0] for x in X])
        ax.set_xticklabels(["$x_{" + str(x[1]) + "," + str(x[2]) + "}$" for x in X])
        ax.tick_params(axis='x', pad=5, labelsize=10)
        ax.tick_params(axis='y', labelsize=10)
        ax.set_xlim(0, 1)
        # plt.yticks([ 0, 0.15])
        # ax.set_xticklabels([0, 1, 2])
        # ax.set_xlabel("index")

        ax.set_yticks([0, 1])

        # copy paste for the 0 centered axis
        # ax.spines['left'].set_position('zero')
        ax.spines['right'].set_color('none')
        ax.spines['bottom'].set_position('zero')
        ax.spines['top'].set_color('none')

        # ax.spines['left'].set_smart_bounds(True)
        # ax.spines['bottom'].set_smart_bounds(True)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_ylabel("$W_" + str(l) + "$", fontsize=10, labelpad=5, rotation=0)

    f.tight_layout()

    plt.show()

    f.savefig("nodale_lineare_basis" + '.pdf', format='pdf', dpi=900, bbox_inches='tight')


def bspline(x,degree):

    if(degree==0):
        if(x>=0 and x<1):
            return 1
        else:
            return 0
    else:

        return (1.0/degree)*(x*bspline(x,degree-1)+(degree+1-x)*bspline(x-1,degree-1))



def bspline_plot(degree):
    fig = plt.figure()

    ax = fig.add_subplot(111)

    font = {'family': 'normal',
            'size': 14}
    plt.rcParams.update({'font.size': 14})

    for p in range(degree+1):



        if(p%2!=0):
            X=np.linspace(0,p+1,500)
            Y=[bspline(x,p) for x in X]
            color = next(ax._get_lines.prop_cycler)['color']
            ax.plot(X, Y,color=color)
            ax.text(p/2+1.3, 0.85-0.05*p, "p="+str(p),color=color)
    ax.spines['bottom'].set_position('zero')
    ax.spines['left'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.set_yticks([0, 1])

    fig.tight_layout()
    fig.savefig("b-splines" + '.pdf', format='pdf', dpi=900, bbox_inches='tight')





X=np.linspace(0,4,500)

Y=[bspline(x,2) for x in X]
print(Y)


#bspline_plot(7)


linearnodalgrid()

#from matplotlib2tikz import save as tikz_save
#tikz_save('nodale_lineare_basis.tex',figureheight='3cm', figurewidth='10cm')


plt.show()