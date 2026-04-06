import logging
import os

def setup_logger(name="indecomposables", log_file=None, verbose=False, debug=False):
    """
    Set up a logger with separate formatters for console and file output.
    
    Args:
        name: Logger name (default: "indecomposables")
        log_file: Path to log file (if None, no file logging)
        verbose: If True, output INFO level messages to stdout
        debug: If True, output DEBUG level messages to log file
    
    Behavior:
        - Console output (stdout):
          * Default: WARNING and above (no timestamp/levelname, message only)
          * With --verbose: INFO and above (message only)
        - File output (if log_file provided):
          * Default: INFO and above (with timestamp/levelname)
          * With --debug: DEBUG and above (with timestamp/levelname)
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    # Avoid duplicate handlers (important in Sage / repeated runs)
    if logger.hasHandlers():
        logger.handlers.clear()

    # Formatter for console: message only (no timestamp/levelname)
    console_formatter = logging.Formatter("%(message)s")

    # Formatter for file: include timestamp and levelname
    file_formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%H:%M:%S"
    )

    # Console handler
    ch = logging.StreamHandler()
    
    # Determine console level based on verbose flag
    if verbose:
        console_level = logging.INFO
    else:
        console_level = logging.WARNING

    ch.setLevel(console_level)
    ch.setFormatter(console_formatter)
    logger.addHandler(ch)

    # File handler
    if log_file:
        # Create parent directories if they don't exist
        log_dir = os.path.dirname(log_file)
        if log_dir and not os.path.exists(log_dir):
            os.makedirs(log_dir)
        
        fh = logging.FileHandler(log_file)
        
        # Determine file level based on debug flag
        if debug:
            file_level = logging.DEBUG
        else:
            file_level = logging.INFO
        
        fh.setLevel(file_level)
        fh.setFormatter(file_formatter)
        logger.addHandler(fh)

    return logger