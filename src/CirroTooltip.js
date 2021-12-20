import {useEffect, useRef} from 'react';

export const SCATTER_TRANSITION = {
    shown: 'opacity 0.2s cubic-bezier(0.23, 1, 0.32, 1) 0s, visibility 0.2s cubic-bezier(0.23, 1, 0.32, 1) 0s, transform 0.4s cubic-bezier(0.23, 1, 0.32, 1) 0s',
    hidden: 'opacity 0.2s cubic-bezier(0.23, 1, 0.32, 1) 0s, transform 0.4s cubic-bezier(0.23, 1, 0.32, 1) 0s'
};


export const DEFAULT_STYLE = {

    position: 'absolute',
    display: 'block',
    borderStyle: 'solid',
    whiteSpace: 'nowrap',
    zIndex: 9999999,
    fontSize: '14px',
    fontFamily: '"Roboto","Helvetica","Arial",sans-serif',
    boxShadow: 'rgba(0, 0, 0, 0.2) 1px 2px 10px',
    backgroundColor: 'rgba(255, 255, 255, 0.9)',
    borderWidth: '1px',
    borderRadius: '4px',
    color: 'rgb(0, 0, 0)',
    padding: '4px',
    top: '0px',
    left: '0px',
    borderColor: 'rgb(84, 112, 198)',
    pointerEvents: 'none',
    visibility: 'hidden'
};


function getTooltipPosition(tip, clientX, clientY, boundsCheck, offset) {
    const parentRect = tip.parentElement.getBoundingClientRect();
    const tipRect = tip.getBoundingClientRect();
    const tipWidth = tipRect.width;
    const tipHeight = tipRect.height;
    const _offset = offset == null ? 10 : offset;
    let left = clientX - parentRect.left + _offset;
    let top = clientY - parentRect.top + _offset;
    // default position is bottom-right

    if (boundsCheck) {
        const scrollBarSize = 18;
        if ((left + tipWidth) >= (parentRect.right - parentRect.left - scrollBarSize)) { // offscreen
            // right, place tip on left
            left = clientX - parentRect.left - _offset - tipWidth;
        }
        if ((top + tipHeight) >= (parentRect.bottom - parentRect.top - scrollBarSize)) { // offscreen
            // bottom, place tip on top
            top = clientY - parentRect.top - _offset - tipHeight;
        }
    }
    return {left: Math.round(left) + 'px', top: Math.round(top) + 'px'};
}

export function setTooltipPosition(tip, left, top, transition) {
    const html = tip.innerHTML;
    tip.style.visibility = html === '' ? 'hidden' : '';
    if (html === '') {
        tip.style.willChange = '';
        if (transition && tip.dataset.shown) {
            tip.style.transition = transition.hidden;
        }
    } else {
        tip.style.willChange = 'transform';
        if (transition) {
            if (!tip.dataset.shown) {
                tip.dataset.shown = 'true';
            } else {
                tip.style.transition = transition.shown;
            }
        }
        tip.style.transform = 'translate3d(' + left + ',' + top + ', 0px)';
    }

}

export default function CirroTooltip(props) {
    const {html, clientX, clientY, offset, style, boundsCheck, transition} = props;
    const _style = Object.assign({}, DEFAULT_STYLE, style);
    const ref = useRef();
    useEffect(() => {
        const tip = ref.current;
        tip.innerHTML = html;
        const position = getTooltipPosition(tip, clientX, clientY, boundsCheck, offset);
        setTooltipPosition(tip, position.left, position.top, transition);
    }, [html, clientX, clientY, boundsCheck, offset, transition]);

    return <div ref={ref} style={_style}></div>;
}

CirroTooltip.defaultProps = {boundsCheck: true};