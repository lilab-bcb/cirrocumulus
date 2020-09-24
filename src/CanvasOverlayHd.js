class CanvasOverlayHd {
    constructor(viewer, options) {
        this._viewer = viewer;
        this.backingScale = 1;
        this._containerWidth = 0;
        this._containerHeight = 0;


        // this._canvasdiv.style.left = 0;
        // this._canvasdiv.style.top = 0;


        this._canvas = document.createElement('canvas');
        this._canvas.style.position = 'absolute';
        this._canvas.style.width = '100%';
        this._canvas.style.height = '100%';
        this._viewer.canvas.appendChild(this._canvas);
        this.onRedraw = options.onRedraw || function () {
        };
        this.clearBeforeRedraw = (typeof (options.clearBeforeRedraw) !== 'undefined') ?
            options.clearBeforeRedraw : true;

        this._viewer.addHandler('update-viewport', () => {
            this.resize();
            this._updateCanvas();
        });
        this._viewer.addHandler('open', () => {
            this.resize();
            this._updateCanvas();
        });
    }

    static getTileIndexFromPixel(viewer, webPoint) {
        let viewportPos = viewer.viewport.pointFromPixel(webPoint);
        for (let i = 0, count = viewer.world.getItemCount(); i < count; i++) {
            let tiledImage = viewer.world.getItemAt(i);
            // box = viewer.world.getItemAt(i).getBounds();
            // if (viewportPos.x > box.x &&
            //   viewportPos.y > box.y &&
            //   viewportPos.x < box.x + box.width &&
            //   viewportPos.y < box.y + box.height) {
            //
            // }
            // tiledImage.lastDrawn.forEach(function (tile) {
            //   if (tile.bounds.containsPoint(viewportPos)) {
            //     console.log('lastDrawn', tile);
            //   }
            // });

            let viewportPosRect = new window.OpenSeadragon.Rect(viewportPos.x, viewportPos.y, 0, 0);
            let tileSourcePosRect = tiledImage._viewportToTiledImageRectangle(viewportPosRect);
            let tileSourcePos = tileSourcePosRect.getTopLeft();
            let source = tiledImage.source;
            if (tileSourcePos.x >= 0 && tileSourcePos.x <= 1 && tileSourcePos.y >= 0 &&
                tileSourcePos.y <= 1 / source.aspectRatio) {
                return i;
                // for (let level = source.minLevel; level <= source.maxLevel; level++) {
                //   let tilePoint = source.getTileAtPoint(level, tileSourcePos);
                //   return i;
                // }
            }
        }
        return -1;
    }

    canvas() {
        return this._canvas;
    }

    context2d() {
        return this._canvas.getContext('2d');
    }

    clear() {
        this._canvas.getContext('2d').clearRect(0, 0, this._containerWidth * this.backingScale, this._containerHeight * this.backingScale);
    }

    resize() {
        let backingScale = 1;
        if (typeof window !== 'undefined' && 'devicePixelRatio' in window) {
            backingScale = window.devicePixelRatio;
        }
        let backingScaleUpdated = this.backingScale !== backingScale;
        this.backingScale = backingScale;
        if (this._containerWidth !== this._viewer.container.clientWidth || backingScaleUpdated) {
            this._containerWidth = this._viewer.container.clientWidth;
            this._canvas.setAttribute('width', backingScale * this._containerWidth);
            this._canvas.setAttribute('width', backingScale * this._containerWidth);
            // this._canvas.style.width = this._containerWidth + 'px';
        }

        if (this._containerHeight !== this._viewer.container.clientHeight || backingScaleUpdated) {
            this._containerHeight = this._viewer.container.clientHeight;
            this._canvas.setAttribute('height', backingScale * this._containerHeight);
            this._canvas.setAttribute('height', backingScale * this._containerHeight);
            // this._canvas.style.height = this._containerHeight + 'px';
        }
    }

    _updateCanvas() {
        let viewportZoom = this._viewer.viewport.getZoom(true);
        if (this.clearBeforeRedraw) {
            this.clear();
        }
        let context = this._canvas.getContext('2d');
        for (let i = 0, count = this._viewer.world.getItemCount(); i < count; i++) {
            let image = this._viewer.world.getItemAt(i);
            if (image) {
                let zoom = image.viewportToImageZoom(viewportZoom);
                let vp = image.imageToViewportCoordinates(0, 0, true);
                let p = this._viewer.viewport.pixelFromPoint(vp, true);
                context.scale(this.backingScale, this.backingScale);
                context.translate(p.x, p.y);
                context.scale(zoom, zoom);
                this.onRedraw({index: i, context: context, x: p.x, y: p.y, zoom: zoom});
                context.setTransform(1, 0, 0, 1, 0, 0);
            }
        }
    }
}

export default CanvasOverlayHd;
