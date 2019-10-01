import pickle as pickle
import os


def createFileName(gridType, model, refineType, maxPoints, maxLevel, degree, objFunc):
    if refineType == 'regular':
        saveName = '{}_{}{}_{}{}.pkl'.format(
            objFunc.getName(), refineType, maxLevel, gridType, degree)
    elif refineType == 'regularByPoints':
        saveName = '{}_{}{}_{}{}.pkl'.format(
            objFunc.getName(), refineType, maxPoints, gridType, degree)
    elif refineType == 'mc':
        saveName = '{}_mc_{}.pkl'.format(objFunc.getName(), maxPoints)
    elif refineType == 'surplus':
        saveName = '{}_{}{}_{}{}.pkl'.format(
            objFunc.getName(), refineType, maxPoints, gridType, degree)
    return saveName


def saveData(data, gridType, model, refineType, maxPoints, maxLevel, degree, objFunc):
    directory = os.path.join(
        '/home/rehmemk/git/SGpp/MR_Python/Extended/data/', model, objFunc.getName())
    if not os.path.exists(directory):
        os.makedirs(directory)

    saveName = createFileName(
        gridType, model, refineType, maxPoints, maxLevel, degree, objFunc)
    datapath = os.path.join(directory, saveName)
    with open(datapath, 'wb') as fp:
        pickle.dump(data, fp)
        print('saved data to {}'.format(datapath))


def loadData(gridType, model, refineType, maxPoints, maxLevel, degree, objFunc):
    directory = os.path.join(
        '/home/rehmemk/git/SGpp/MR_Python/Extended/data/', model, objFunc.getName())
    saveName = createFileName(
        gridType, model, refineType, maxPoints, maxLevel, degree, objFunc)
    datapath = os.path.join(directory, saveName)
    with open(datapath, 'rb') as fp:
        # encoding necessary because the data was stored with cPickle from python2.7
        # data = pickle.load(fp, encoding='latin1')
        data = pickle.load(fp, encoding='latin1')
    print('loaded data from {}'.format(datapath))
    return data
