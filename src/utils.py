#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Imports
import json
import numpy as np
import pandas as pd
import os


def load_data_from_json(filename):
    """
    Load data from a JSON file.

    Parameters:
        filename (str): Path to the JSON file.

    Returns:
        dict: Data loaded from the JSON file.

    Raises:
        FileNotFoundError: If the file does not exist.
        json.JSONDecodeError: If the file contains invalid JSON.
    """
    try:
        with open(filename, 'r') as file:
            return json.load(file)
    except json.JSONDecodeError as e:
        raise json.JSONDecodeError(
            f"Error decoding JSON from {filename}: {e.msg}",
            e.doc, e.pos
        )
    except FileNotFoundError as e:
        raise FileNotFoundError(f"File not found: {filename}") from e


def save_data_as_json(data, filename, indent=4):
    """
    Save data to a JSON file, converting non-serializable keys and values to
    appropriate JSON-compatible types.

    Parameters:
        data: The data to save.
        filename (str): The filename for the JSON file.
        indent (int, optional): Number of spaces for indentation. Default is 4.

    Raises:
        OSError: If there's an error writing to the file.
        TypeError: If data contains types that cannot be converted.
    """
    try:
        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(os.path.abspath(filename)), exist_ok=True)
        
        # Convert data to JSON-serializable format
        data_serializable = _convert_for_json(data)
        
        # Write to file
        with open(filename, 'w') as file:
            json.dump(data_serializable, file, indent=indent)

    except (OSError, TypeError) as e:
        raise type(e)(f"Error saving data to {filename}: {str(e)}") from e


def _is_numeric_key_dict(obj):
    """Check if dict has numeric string keys (like positions)."""
    if not isinstance(obj, dict) or len(obj) == 0:
        return False
    try:
        # Check if all keys can be converted to int
        for key in obj.keys():
            int(str(key))
        return True
    except (ValueError, TypeError):
        return False


def _convert_for_json(obj):
    """
    Convert Python objects to JSON-serializable types.

    Parameters:
        obj: Any Python object.

    Returns:
        Object converted to JSON-serializable type.
    """
    if isinstance(obj, dict):
        # Sort by numeric key if dict contains position-like keys
        if _is_numeric_key_dict(obj):
            sorted_items = sorted(obj.items(), key=lambda x: int(str(x[0])))
            return {str(key): _convert_for_json(value) for key, value in sorted_items}
        else:
            return {str(key): _convert_for_json(value) for key, value in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [_convert_for_json(element) for element in obj]
    elif isinstance(obj, set):
        return list(obj)  # Convert sets to lists
    elif isinstance(obj, np.ndarray):
        return obj.tolist()  # Convert NumPy arrays to lists
    elif isinstance(obj, (np.integer, np.floating)):
        return obj.item()  # Convert NumPy numbers to Python native types
    elif isinstance(obj, np.bool_):
        return bool(obj)  # Convert NumPy booleans to Python bool
    elif pd.isna(obj):  # Handle pandas NA/NaN values
        return None
    elif isinstance(obj, pd.DataFrame):
        return obj.to_dict('records')  # Convert DataFrame to list of records
    elif isinstance(obj, pd.Series):
        return obj.to_dict()  # Convert Series to dict
    else:
        return obj  # Return unchanged if already serializable