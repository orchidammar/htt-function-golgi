"""
Stage 2: Cell Segmentation

Purpose:
    Split each image into individual cells using C1 as reference.

CRITICAL: C1-as-Reference Strategy
    - Segment C1 (before/baseline) to detect N cells
    - Use those N cell bounding boxes as reference
    - Extract same spatial regions from C2 and C3
    - This handles the fact that particles scatter in C2/C3

Input:
    - outputs/01_validated/manifest.json
    - Raw images from test_data/
    - Config: segmentation parameters

Output:
    - outputs/02_segmented/[group]/[sample_id]/C1/cell_[n].tif
    - outputs/02_segmented/[group]/[sample_id]/C2/cell_[n].tif
    - outputs/02_segmented/[group]/[sample_id]/C3/cell_[n].tif
    - outputs/02_segmented/[group]/[sample_id]/sample_metadata.json
"""

import json
import cv2
import numpy as np
from pathlib import Path
from datetime import datetime
from tqdm import tqdm
import click


class CellSegmenter:
    """Segments cells using C1-as-reference strategy"""

    def __init__(self, config, manifest, output_dir):
        self.config = config
        self.manifest = manifest
        self.output_dir = Path(output_dir)

        # Segmentation parameters from config
        seg_config = config.get('segmentation', {})
        self.min_cell_area = seg_config.get('min_cell_area', 100)
        self.max_cell_area = seg_config.get('max_cell_area', 50000)
        self.threshold_method = seg_config.get('threshold_method', 'otsu')

        # Statistics
        self.stats = {
            'samples_processed': 0,
            'samples_failed': 0,
            'total_cells_detected': 0,
            'cells_excluded': 0,
            'cells_valid': 0
        }

        self.segmentation_report = {
            'created_at': datetime.now().isoformat(),
            'samples': {},
            'errors': [],
            'warnings': []
        }

    def run(self):
        """Execute Stage 2: Cell Segmentation"""
        click.echo("\n" + "=" * 60)
        click.echo("[Stage 2] Cell Segmentation")
        click.echo("=" * 60)

        samples = self.manifest['samples']
        valid_samples = {sid: s for sid, s in samples.items() if s.get('valid', False)}

        click.echo(f"\n→ Segmenting {len(valid_samples)} samples using C1-as-reference strategy...")

        for sample_id, sample_data in tqdm(valid_samples.items(), desc="Segmenting"):
            try:
                self._segment_sample(sample_id, sample_data)
                self.stats['samples_processed'] += 1
            except Exception as e:
                click.echo(f"\n✗ Error segmenting {sample_id}: {str(e)}")
                self.stats['samples_failed'] += 1
                self.segmentation_report['errors'].append({
                    'sample': sample_id,
                    'error': str(e)
                })

        # Save segmentation report
        self._save_report()

        # Print results
        self._print_results()

        return self.segmentation_report

    def _segment_sample(self, sample_id, sample_data):
        """Segment a single sample using C1-as-reference strategy"""

        # Create output directories
        group = sample_data['group']
        sample_name = sample_data['sample_name']
        sample_output_dir = self.output_dir / group / sample_name

        for channel in sample_data['channels'].keys():
            (sample_output_dir / channel).mkdir(parents=True, exist_ok=True)

        # Step 1: Load C1 image (reference channel)
        c1_path = sample_data['channels']['C1']['path']
        c1_img = cv2.imread(c1_path, cv2.IMREAD_UNCHANGED)

        if c1_img is None:
            raise ValueError(f"Failed to load C1 image: {c1_path}")

        # Step 2: Detect cells in C1
        cells = self._detect_cells_in_c1(c1_img, sample_id)

        if len(cells) == 0:
            self.segmentation_report['warnings'].append({
                'code': 'E302',
                'sample': sample_id,
                'message': 'No cells detected in C1'
            })
            return

        # Step 3: For each detected cell, extract from all channels
        valid_cells = []

        for cell_idx, cell_info in enumerate(cells):
            bbox = cell_info['bbox']
            x, y, w, h = bbox

            # Validate cell (check boundaries, size, etc.)
            if not self._validate_cell(cell_info, c1_img.shape, sample_id, cell_idx):
                self.stats['cells_excluded'] += 1
                continue

            # Extract cell from all channels using SAME bbox
            cell_id = f"cell_{cell_idx + 1}"
            cell_metadata = {
                'cell_id': cell_id,
                'bbox': bbox,
                'center': cell_info['center'],
                'area': cell_info['area'],
                'channels': {}
            }

            for channel, channel_data in sample_data['channels'].items():
                # Load channel image
                img = cv2.imread(channel_data['path'], cv2.IMREAD_UNCHANGED)

                if img is None:
                    raise ValueError(f"Failed to load {channel} image")

                # Extract cell using same bbox from C1
                cell_img = img[y:y+h, x:x+w].copy()

                # Save cell image
                cell_path = sample_output_dir / channel / f"{cell_id}.tif"
                cv2.imwrite(str(cell_path), cell_img)

                # Save metadata
                cell_metadata['channels'][channel] = {
                    'path': str(cell_path),
                    'shape': cell_img.shape,
                    'mean_intensity': float(np.mean(cell_img)),
                    'max_intensity': float(np.max(cell_img)),
                    'min_intensity': float(np.min(cell_img))
                }

            valid_cells.append(cell_metadata)
            self.stats['cells_valid'] += 1

        # Save sample metadata
        sample_metadata = {
            'sample_id': sample_id,
            'group': group,
            'sample_name': sample_name,
            'cells_detected': len(cells),
            'cells_valid': len(valid_cells),
            'cells_excluded': len(cells) - len(valid_cells),
            'cells': valid_cells,
            'segmentation_params': {
                'min_area': self.min_cell_area,
                'max_area': self.max_cell_area,
                'threshold_method': self.threshold_method
            }
        }

        metadata_path = sample_output_dir / 'sample_metadata.json'
        with open(metadata_path, 'w') as f:
            json.dump(sample_metadata, f, indent=2)

        self.segmentation_report['samples'][sample_id] = sample_metadata

    def _detect_cells_in_c1(self, img, sample_id):
        """Detect cells in C1 using contour-based detection"""

        # Convert to grayscale if RGB
        if len(img.shape) == 3:
            gray = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)
        else:
            gray = img.copy()

        # Preprocessing
        # 1. Gaussian blur to reduce noise
        blurred = cv2.GaussianBlur(gray, (5, 5), 0)

        # 2. Morphological operations
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))

        # Close (fill small holes)
        closed = cv2.morphologyEx(blurred, cv2.MORPH_CLOSE, kernel, iterations=2)

        # Open (remove small noise)
        opened = cv2.morphologyEx(closed, cv2.MORPH_OPEN, kernel, iterations=1)

        # 3. Thresholding
        if self.threshold_method == 'otsu':
            threshold_value, binary = cv2.threshold(
                opened, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU
            )
        else:
            # Default to Otsu if method unknown
            threshold_value, binary = cv2.threshold(
                opened, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU
            )

        # 4. Find contours
        contours, _ = cv2.findContours(
            binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE
        )

        # 5. Filter contours by area and extract cell info
        cells = []

        for contour in contours:
            area = cv2.contourArea(contour)

            # Filter by area
            if area < self.min_cell_area or area > self.max_cell_area:
                continue

            # Get bounding box
            x, y, w, h = cv2.boundingRect(contour)

            # Calculate center
            M = cv2.moments(contour)
            if M['m00'] != 0:
                cx = int(M['m10'] / M['m00'])
                cy = int(M['m01'] / M['m00'])
            else:
                cx = x + w // 2
                cy = y + h // 2

            # Calculate shape metrics
            perimeter = cv2.arcLength(contour, True)
            circularity = 4 * np.pi * area / (perimeter * perimeter) if perimeter > 0 else 0
            aspect_ratio = max(w, h) / min(w, h) if min(w, h) > 0 else 0

            cells.append({
                'bbox': [x, y, w, h],
                'center': [cx, cy],
                'area': area,
                'perimeter': perimeter,
                'circularity': circularity,
                'aspect_ratio': aspect_ratio,
                'contour': contour
            })

        # Sort by area (largest first)
        cells.sort(key=lambda c: c['area'], reverse=True)

        self.stats['total_cells_detected'] += len(cells)

        return cells

    def _validate_cell(self, cell_info, img_shape, sample_id, cell_idx):
        """Validate a detected cell

        Excludes cells within 10% of image edge (boundary exclusion for PhD thesis rigor)
        """

        x, y, w, h = cell_info['bbox']
        img_h, img_w = img_shape[:2]

        # Get cell center
        cx, cy = cell_info['center']

        # Boundary exclusion: 10% of image edge
        # Get from config, default to 0.1 (10%)
        boundary_threshold = self.config.get('validation', {}).get('boundary_exclusion_threshold', 0.1)
        margin_x = int(img_w * boundary_threshold)
        margin_y = int(img_h * boundary_threshold)

        # Check if cell center is within boundary zone
        if (cx < margin_x or cx > (img_w - margin_x) or
            cy < margin_y or cy > (img_h - margin_y)):
            self.segmentation_report['warnings'].append({
                'code': 'E303',
                'sample': sample_id,
                'cell': cell_idx,
                'message': f'Cell at boundary (center within {boundary_threshold*100:.0f}% edge zone), excluded',
                'boundary_threshold': boundary_threshold,
                'cell_center': [cx, cy],
                'image_margins': [margin_x, margin_y, img_w - margin_x, img_h - margin_y]
            })
            return False

        # Check cell size
        area = cell_info['area']
        if area < self.min_cell_area:
            self.segmentation_report['warnings'].append({
                'code': 'E305',
                'sample': sample_id,
                'cell': cell_idx,
                'message': f'Cell too small ({area} px²), excluded'
            })
            return False

        if area > self.max_cell_area:
            self.segmentation_report['warnings'].append({
                'code': 'E306',
                'sample': sample_id,
                'cell': cell_idx,
                'message': f'Cell too large ({area} px²), flagged but kept'
            })
            # Don't exclude, just flag

        # Check aspect ratio (elongated cells)
        aspect_ratio = cell_info['aspect_ratio']
        if aspect_ratio > 8:
            self.segmentation_report['warnings'].append({
                'code': 'E308',
                'sample': sample_id,
                'cell': cell_idx,
                'message': f'Invalid cell shape (aspect ratio {aspect_ratio:.1f}), flagged'
            })
            # Don't exclude, just flag

        return True

    def _save_report(self):
        """Save segmentation report"""
        report_path = self.output_dir / 'segmentation_report.json'

        self.segmentation_report['summary'] = {
            'samples_processed': self.stats['samples_processed'],
            'samples_failed': self.stats['samples_failed'],
            'total_cells_detected': self.stats['total_cells_detected'],
            'cells_valid': self.stats['cells_valid'],
            'cells_excluded': self.stats['cells_excluded'],
            'exclusion_rate': (self.stats['cells_excluded'] / self.stats['total_cells_detected'] * 100)
                if self.stats['total_cells_detected'] > 0 else 0
        }

        with open(report_path, 'w') as f:
            json.dump(self.segmentation_report, f, indent=2)

    def _print_results(self):
        """Print segmentation results"""
        summary = self.segmentation_report['summary']

        click.echo(f"\n✓ Stage 2 completed")
        click.echo(f"  → Output: {self.output_dir}")
        click.echo(f"\n  Samples processed: {summary['samples_processed']}")
        click.echo(f"  Samples failed: {summary['samples_failed']}")
        click.echo(f"  Total cells detected: {summary['total_cells_detected']}")
        click.echo(f"  Cells valid: {summary['cells_valid']}")
        click.echo(f"  Cells excluded: {summary['cells_excluded']} ({summary['exclusion_rate']:.1f}%)")

        # Count boundary exclusions specifically
        boundary_exclusions = sum(1 for w in self.segmentation_report['warnings'] if w.get('code') == 'E303')
        if boundary_exclusions > 0:
            click.echo(f"    → Boundary exclusions: {boundary_exclusions} (within 10% of image edge)")

        if len(self.segmentation_report['errors']) > 0:
            click.echo(f"\n⚠ {len(self.segmentation_report['errors'])} errors. Check segmentation_report.json")

        if len(self.segmentation_report['warnings']) > 0:
            click.echo(f"⚠ {len(self.segmentation_report['warnings'])} warnings. Check segmentation_report.json")


def run_stage_2(config, manifest, output_dir):
    """Run Stage 2: Cell Segmentation"""
    segmenter = CellSegmenter(config, manifest, output_dir)
    return segmenter.run()
