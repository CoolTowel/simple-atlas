import time
import multiprocessing
from functools import wraps

from skychart import *

def timefunc(func):
    '''
    time the given function
    '''
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        func_result = func(*args, **kwargs)
        end_time = time.time()
        print("function {} takes {}".format(func.__name__, end_time-start_time))
        return func_result
    return wrapper

def skymap_save(fig, filename, echo):
    if echo:
        print("Saving to "+filename)
    fig.savefig(filename)
    plt.close(fig)


@timefunc
def plot_savefig_synchronous(num_fig=10, echo=False):
    for i, t in enumerate(time_ini + np.arange(num_fig) * u.day):
        fig = skychart(time=t, location=location, show=False)
        save_path = './test/ani_365/{}.png'.format(i)
        skymap_save(fig, save_path, echo)

@timefunc
def plot_savefig_async(num_fig=10, echo=False):
    processes = []
    for i, t in enumerate(time_ini + np.arange(num_fig) * u.day):
        fig = skychart(time=t, location=location, show=False)
        save_path = './test/ani_365/{}.png'.format(i)
        # Save figure
        p = multiprocessing.Process(target=skymap_save, args=(fig, save_path, echo,))
        p.start()
        processes.append(p)
    for p in processes:
        p.join()



if __name__ == '__main__':
    time_ini = Time('2022-1-1 21:00:00')
    loc_lat = 52.15485022021291
    loc_lon = 4.483882499710578
    loc_height = 0
    location = EarthLocation(lat=loc_lat * u.deg,
                             lon=loc_lon * u.deg,
                             height=loc_height * u.m)

    # Benchmarking
    test_size = 300
    echo = False
    plot_savefig_synchronous(num_fig=test_size, echo=echo)
    plot_savefig_async(num_fig=test_size, echo=echo)