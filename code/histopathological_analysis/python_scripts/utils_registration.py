import os
import pathlib
from valis import registration, warp_tools, viz, slide_tools

class customValis(registration.Valis):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def create_non_rigid_reg_mask(self):
        """
        Get mask for non-rigid registration
        """
        non_rigid_mask = self.mask_dict["all"][0]

        for slide_obj in self.slide_dict.values():
            slide_obj.non_rigid_reg_mask = non_rigid_mask

        # Save thumbnail of mask
        ref_slide = self.get_ref_slide()
        if ref_slide.img_type == slide_tools.IHC_NAME:
            ref_img = warp_tools.resize_img(ref_slide.image, ref_slide.processed_img_shape_rc)
            warped_ref_img = ref_slide.warp_img(ref_img, non_rigid=False, crop="reference")
        else:
            warped_ref_img = ref_slide.warp_img(ref_slide.processed_img, non_rigid=False, crop="reference")

        pathlib.Path(self.mask_dir).mkdir(exist_ok=True, parents=True)
        thumbnail_img = self.create_thumbnail(warped_ref_img)

        draw_mask = warp_tools.resize_img(non_rigid_mask, ref_slide.reg_img_shape_rc, interp_method="nearest")
        _, overlap_mask_bbox_xywh = self.get_crop_mask("reference")
        draw_mask = warp_tools.crop_img(draw_mask, overlap_mask_bbox_xywh.astype(int))
        thumbnail_mask = self.create_thumbnail(draw_mask)

        thumbnail_mask_outline = viz.draw_outline(thumbnail_img, thumbnail_mask)
        outline_f_out = os.path.join(self.mask_dir, f'{self.name}_non_rigid_mask.png')
        warp_tools.save_img(outline_f_out, thumbnail_mask_outline)