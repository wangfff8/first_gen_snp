import logging
import sys


class ColoredFormatter(logging.Formatter):
    """Custom logging formatter to provide colored output based on level."""

    grey = "\x1b[38;20m"
    green = "\x1b[32;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"

    # Define different formats for different levels
    FORMATS = {
        logging.DEBUG: grey + "%(levelname)s: %(message)s" + reset,
        logging.INFO: green + "%(message)s" + reset,
        logging.WARNING: yellow + "WARNING: %(message)s" + reset,
        logging.ERROR: red + "ERROR: %(message)s" + reset,
        logging.CRITICAL: bold_red + "CRITICAL: %(message)s" + reset,
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


def setup_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    """Sets up a logger with colored console output if connected to a TTY."""
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # Avoid adding multiple handlers if the logger already has some.
    if not logger.handlers:
        handler = logging.StreamHandler()
        if sys.stdout.isatty():
            handler.setFormatter(ColoredFormatter())
        else:
            handler.setFormatter(logging.Formatter("%(message)s"))

        logger.addHandler(handler)
        logger.propagate = False

    return logger
