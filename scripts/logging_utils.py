import logging

def setup_logger(name="indecomposables", log_file=None, verbose=False, debug=False):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    # Avoid duplicate handlers (important in Sage / repeated runs)
    if logger.hasHandlers():
        logger.handlers.clear()

    # Formatter
    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%H:%M:%S"
    )

    # Console handler
    ch = logging.StreamHandler()
    if debug:
        console_level = logging.DEBUG
    elif verbose:
        console_level = logging.INFO
    else:
        console_level = logging.WARNING

    ch.setLevel(console_level)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # File handler
    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger