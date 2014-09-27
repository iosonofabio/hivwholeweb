#!/usr/bin/python2

def static_plot_example():
    #import sys
    #sys.path.append('/usr/lib/python2.7/site-packages')

    import numpy as np
    from bokeh.plotting import line
    import bokeh.embed as embed
    from bokeh.resources import CDN

    x_data = np.arange(1, 101)
    y_data = np.random.randint(0, 101, 100)

    plot = line(x_data, y_data)

    js, div = embed.components(plot, CDN)

    return js, div


def coverage(fragment='F1'):
    import sys
    hivwholeseq_path = '/home/fabio/university/phd/sequencing/script/mapping'
    if hivwholeseq_path not in sys.path:
        sys.path.append(hivwholeseq_path)
    general_path = '/usr/lib/python2.7/site-packages'
    if general_path not in sys.path:
        sys.path.append(general_path)
    
    from hivwholeseq.patients.patients import load_sample_sequenced
    sample = load_sample_sequenced('VL96-15555')
    cov = sample.get_coverage(fragment=fragment, PCR=1)

    import numpy as np
    from bokeh.plotting import line
    import bokeh.embed as embed
    from bokeh.resources import CDN

    x_data = np.arange(cov.shape[-1])
    y_data = cov
    plot = line(x_data, y_data, title=sample.name+', '+fragment)

    js, div = embed.components(plot, CDN)

    return js, div


def build_dynamic_plot():
    import numpy as np
    #from bokeh.plotting import *
    import bokeh.embed as embed
    from bokeh.resources import CDN

    return js, tag


#def parse_prepare_bokeh_stdout(stream):
#    '''Parse the output of the stdout stream'''
#    tag = ''
#    for line in iter(stream.readline, ''):
#        tag += line
#        print tag
#        if '</script>' in line:
#            break
#
#    print tag
#
#    tag = tag.strip('\n')
#    taglines = tag.replace(' ', '').split('\n')
#    for i, line in enumerate(taglines):
#        if 'START_JAVASCRIPT' in line:
#            start_js = i
#        elif 'END_JAVASCRIPT' in line:
#            end_js = i
#        elif '<script' in line:
#            taglines = taglines[i + 1: - 1]
#            break
#
#    js = taglines[start_js:end_js]
#    print '\n'.join(js)
#
#    context_dict = {'bokeh_js': js}
#    for f in taglines:
#        print f
#        key, value = f.split('=')
#        key = key.replace('-', '_')
#        value = value.strip('"')
#        context_dict[key] = value
#
#    return context_dict


#def spawn_bokeh_client(pop, timeout=50, session=None):
#    import time
#    
#    t0 = time.time()
#    while time.time() - t0 < timeout:
#        pop.session.load_document(pop.document)
#        time.sleep(0.5)





if __name__ == '__main__':

    import sys
    args = sys.argv
    if len(args) > 1:
        timeout = int(args[1])
    else:
        timeout = 50
    if len(args) > 2:
        js_static = args[2]
    else:
        js_static = None

    if js_static is None:
        tag, pop = prepare_bokeh_client()
    else:
        js, tag, pop = prepare_bokeh_client(js_static=js_static)
        sys.stdout.write('START_JAVASCRIPT\n')
        sys.stdout.write(js)
        sys.stdout.write('\nEND_JAVASCRIPT\n')

    sys.stdout.write(tag)
    sys.stdout.flush()
    sys.stdout.write('\n')
    sys.stdout.flush()

    #spawn_bokeh_client(pop, timeout=timeout)
