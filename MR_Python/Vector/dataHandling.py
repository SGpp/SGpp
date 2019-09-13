import pickle
import os


def createFileName(gridType, model, refineType, maxPoints, maxLevel, degree, objFunc):
    if refineType == 'regular':
        saveName = '{}_{}{}_{}{}.pkl'.format(
            objFunc.getName(), refineType, maxLevel, gridType, degree)
    elif refineType == 'mc':
        saveName = '{}_mc_{}.pkl'.format(objFunc.getName(), maxPoints)
    elif refineType == 'surplus':
        saveName = '{}_{}{}_{}{}.pkl'.format(
            objFunc.getName(), refineType, maxPoints, gridType, degree)
    return saveName


def saveData(data, path, gridType, model, refineType, maxPoints, maxLevel, degree, objFunc):
    directory = os.path.join(path, model, objFunc.getName())
    if not os.path.exists(directory):
        os.makedirs(directory)

    saveName = createFileName(
        gridType, model, refineType, maxPoints, maxLevel, degree, objFunc)
    datapath = os.path.join(directory, saveName)
    with open(datapath, 'wb') as fp:
        pickle.dump(data, fp)
        print('saved data to {}'.format(datapath))


def loadData(path, gridType, model, refineType, maxPoints, maxLevel, degree, objFunc):
    directory = os.path.join(path, model, objFunc.getName())
    saveName = createFileName(
        gridType, model, refineType, maxPoints, maxLevel, degree, objFunc)
    datapath = os.path.join(directory, saveName)
    with open(datapath, 'rb') as fp:
        data = pickle.load(fp)
    print('loaded data from {}'.format(datapath))
    return data
