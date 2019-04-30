import logging.config

from config import LOG_FILE


def setup_logging(fname=LOG_FILE):
    NO_COLOR = "\33[m"
    RED, GREEN, ORANGE, BLUE, PURPLE, LBLUE, GREY = (
        map("\33[%dm".__mod__, range(31, 38)))

    logging_config = {
        'version':  1,
        'disable_existing_loggers': False,
        'formatters': {
            'default': {
                'format': '[%(asctime)s] %(levelname)-4s - '
                          '%(name)-4s - %(message)s'
            },
            'color': {
                'format': '{}[%(asctime)s]{} {}%(levelname)-5s{} - '
                          '{}%(name)-5s{}: %(message)s'.format(
                              GREEN, NO_COLOR, PURPLE, NO_COLOR,
                              ORANGE, NO_COLOR)
            }
        },
        'handlers': {
            'stream': {
                'class': 'logging.StreamHandler',
                'formatter': 'color',
            }
        },
        'root': {
            'handlers': ['stream'],
            'level': logging.INFO,
        },
    }
    if fname is not None:
        logging_config['handlers']['file'] = {
            'class': 'logging.FileHandler',
            'formatter': 'default',
            'level': logging.DEBUG,
            'filename': fname,
        }
        logging_config['root']['handlers'].append('file')

    logging.config.dictConfig(logging_config)
