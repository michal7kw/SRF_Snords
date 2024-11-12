import os
import sys
import nbformat
import logging
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def setup_scripts_directory():
    """Create Scripts directory if it doesn't exist."""
    scripts_dir = Path('./Scripts')
    try:
        scripts_dir.mkdir(exist_ok=True)
        logger.info(f"Scripts directory ready at {scripts_dir.absolute()}")
        return scripts_dir
    except Exception as e:
        logger.error(f"Failed to create Scripts directory: {e}")
        sys.exit(1)

def convert_notebook_to_python(notebook_path, output_dir):
    """Convert a single .ipynb file to .py format."""
    try:
        # Read the notebook
        with open(notebook_path, 'r', encoding='utf-8') as f:
            notebook = nbformat.read(f, as_version=4)

        # Create output path
        output_path = output_dir / f"{notebook_path.stem}.py"
        
        # Convert and write to .py file
        python_code = []
        
        # Add a header comment
        python_code.append(f"# Converted from {notebook_path.name}\n\n")
        
        for cell in notebook['cells']:
            if cell['cell_type'] == 'code':
                # Add cell source code
                python_code.append(cell['source'])
                # Add newline between cells
                python_code.append('\n\n')
            elif cell['cell_type'] == 'markdown':
                # Convert markdown cells to comments
                commented_text = '\n'.join(f"# {line}" if line.strip() else "#" 
                                         for line in cell['source'].split('\n'))
                python_code.append(commented_text)
                python_code.append('\n\n')

        # Write the Python file
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(''.join(python_code))
            
        logger.info(f"Successfully converted {notebook_path.name} to {output_path.name}")
        return True
        
    except Exception as e:
        logger.error(f"Error converting {notebook_path}: {e}")
        return False

def main():
    """Main function to convert all notebooks in current directory."""
    # Setup Scripts directory
    scripts_dir = setup_scripts_directory()
    
    # Get all .ipynb files in current directory
    notebooks = list(Path('.').glob('*.ipynb'))
    
    if not notebooks:
        logger.warning("No .ipynb files found in current directory")
        return
    
    logger.info(f"Found {len(notebooks)} notebook(s) to convert")
    
    # Convert each notebook
    successful = 0
    failed = 0
    
    for notebook_path in notebooks:
        if convert_notebook_to_python(notebook_path, scripts_dir):
            successful += 1
        else:
            failed += 1
    
    # Print summary
    logger.info(f"\nConversion complete:")
    logger.info(f"Successfully converted: {successful}")
    logger.info(f"Failed conversions: {failed}")

if __name__ == "__main__":
    main()