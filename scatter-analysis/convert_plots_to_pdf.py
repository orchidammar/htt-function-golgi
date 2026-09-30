#!/usr/bin/env python3
"""
Convert all PNG plots to PDF format in separate directory with _pdf suffix.
"""

from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import shutil

def convert_png_to_pdf(png_path, base_source_dir, base_pdf_dir):
    """Convert a single PNG file to PDF in mirrored directory structure"""

    # Calculate relative path from source base dir
    rel_path = png_path.relative_to(base_source_dir)

    # Create corresponding PDF path
    pdf_rel_path = rel_path.with_suffix('.pdf')
    pdf_path = base_pdf_dir / pdf_rel_path

    # Create parent directory if needed
    pdf_path.parent.mkdir(parents=True, exist_ok=True)

    # Read the PNG image
    img = mpimg.imread(png_path)

    # Create figure with exact image dimensions
    dpi = 300
    height, width = img.shape[:2]
    figsize = width / dpi, height / dpi

    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    ax.imshow(img)

    # Save as PDF
    plt.savefig(pdf_path, dpi=dpi, bbox_inches='tight', pad_inches=0)
    plt.close(fig)

    return pdf_path

def main():
    """Find all PNG files in pooled_plots and convert to PDF in pooled_plots_pdf"""
    base_source_dir = Path("results/plots/png")
    base_pdf_dir = Path("results/plots/pdf")

    if not base_source_dir.exists():
        print(f"Error: {base_source_dir} does not exist")
        return

    # Create PDF base directory
    base_pdf_dir.mkdir(exist_ok=True)

    # Find all PNG files recursively
    png_files = list(base_source_dir.rglob("*.png"))

    if not png_files:
        print(f"No PNG files found in {base_source_dir}")
        return

    print(f"\nConverting {len(png_files)} PNG files to PDF")
    print(f"Source: {base_source_dir}/")
    print(f"Target: {base_pdf_dir}/")
    print("="*80)

    converted = []
    for png_path in sorted(png_files):
        rel_path = png_path.relative_to(base_source_dir)
        print(f"Converting: {rel_path}")
        pdf_path = convert_png_to_pdf(png_path, base_source_dir, base_pdf_dir)
        print(f"  → Saved: {pdf_path}")
        converted.append(pdf_path)

    print("\n" + "="*80)
    print(f"✅ Successfully converted {len(converted)} files to PDF")
    print(f"\nAll PDFs saved to: {base_pdf_dir}/")
    print("\nDirectory structure:")

    # Show directory structure
    for pdf_path in sorted(converted):
        print(f"  • {pdf_path}")
    print()

if __name__ == '__main__':
    main()
