import numpy as np


def addMotion(data):
    """This function randomly adds 'motion' to data collected from a
    stationary gyroscope, making it seem as if the data was collected
    from a moving vehicle
    """
    num = len(data)
    count = 15

    dur = np.random.normal(3.5*45, 4*45, (count,))
    w = np.random.normal(0.05, 0.5, (count,))

    dur = np.hstack((dur, np.random.normal(2*45, 3.5*45, (count,))))
    w = np.hstack((w, np.random.normal(1.5, 2, (count,))))

    dur = np.hstack((dur, np.random.normal(1*45, 2*45, (count,))))
    w = np.hstack((w, np.random.normal(3, 4, (count,))))

    idx = np.arange(num)
    np.random.shuffle(idx)

    idx = idx[:len(dur)]

    for i in xrange(len(idx)):
        data[idx[i]:idx[i]+dur[i]] = data[idx[i]:idx[i]+dur[i]] + w[i]

    return data
