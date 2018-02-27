import multiprocessing,os,warnings,sys
from multiprocessing import Pool
import numpy as np

if sys.version_info[0] < 3:
    import cPickle as pick
else:
    import pickle as pick

def foreach(toPar,parFunc,args,numThreads=multiprocessing.cpu_count()):
    p = Pool(processes=numThreads)
    results=[]
    for x in p.imap_unordered(_parWrap,[[parFunc,np.append(y,args)] for y in toPar]):
        results.append(x)
    p.close()

    outddict=dict([])
    if isinstance(results[0],dict):
        for res in results:
            outddict[res['key']]={k:res[k] for k in res.keys() if k != 'key'}
        return outddict
    else:
        return results

def _parWrap(args):
    func,newArgs=args
    try:
        return(func(newArgs))
    except RuntimeError:
        print('something')
        return(None)

def _pickleable(obj):
    try:
        with open(r"temp.pickle", "wb") as output_file:
            pick.dump(obj, output_file)
        pickle=True
    except:
        pickle=False
    try:
        os.remove('temp.pickle')
    except:
        pass
    return pickle

def parReturn(toReturn):
    if isinstance(toReturn,dict):
        final=dict([])
        for key in toReturn:
            if _pickleable(toReturn[key]):
                final[key]=toReturn[key]
            else:
                print("Had to remove object %s from return dictionary, as it was not pickleable."%key)

    elif isinstance(toReturn,(tuple,list,np.array)):
        final=[]
        for i in range(len(toReturn)):
            if _pickleable(toReturn[i]):
                final.append(toReturn[i])
            else:
                warnings.warn(RuntimeWarning,"Had to remove the %i (th) object from return array, as it was not pickleable."%i)
    else:
        print('I do not recognize the data type of your return variable')
        sys.exit()

    if final:
        return final
    else:
        print('Nothing you wanted to return was Pickleable')
        sys.exit()
'''
def parReturn(toReturn):
    name=False
    if isinstance(toReturn,dict):
        final=dict([])
        for key in toReturn:
            if key is 'key':
                name=True
            if _pickleable(toReturn[key]):
                final[key]=toReturn[key]
            else:
                print("Had to remove object %s from return dictionary, as it was not pickleable."%key)

        if not name:
            print('You must have a "key" element of your return dictionary, so that this result dictionary can be identified later.')
            sys.exit()
    elif isinstance(toReturn,(tuple,list,np.array)):
        final=[]
        if not isinstance(toReturn[0],str):
            print('The first element of your return array must be an identifying string.')
            sys.exit()
        for i in range(len(toReturn)):
            if _pickleable(toReturn[i]):
                final.append(toReturn[i])
            else:
                warnings.warn(RuntimeWarning,"Had to remove the %i (th) object from return array, as it was not pickleable."%i)
    else:
        print('I do not recognize the data type of your return variable')
        sys.exit()


    return final
'''